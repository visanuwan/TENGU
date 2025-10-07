import numpy as np
import pandas as pd
import scipy.sparse
import scipy.spatial
import scanpy as sc
import anndata as ad
from copy import deepcopy

def destripe_counts(adata, counts_key="n_counts", adjusted_counts_key="n_counts_adjusted"):
    """Scale each bin's counts to a new total."""
    sc.get._utils.view_to_actual(adata)
    scaling_factors = adata.obs[adjusted_counts_key] / adata.obs[counts_key]
    adata.X = scipy.sparse.diags(scaling_factors).dot(adata.X)

def destripe(adata, quantile=0.99, counts_key="n_counts", factor_key="destripe_factor", adjusted_counts_key="n_counts_adjusted", adjust_counts=True):
    """Correct raw counts for variable bin widths in VisiumHD."""
    quant_row = adata.obs.groupby("array_row")[counts_key].quantile(quantile)
    adata.obs[factor_key] = adata.obs[counts_key] / adata.obs["array_row"].map(quant_row)
    
    quant_col = adata.obs.groupby("array_col")[factor_key].quantile(quantile)
    adata.obs[factor_key] /= adata.obs["array_col"].map(quant_col)
    
    adata.obs[adjusted_counts_key] = adata.obs[factor_key] * np.quantile(adata.obs[counts_key], quantile)
    
    if adjust_counts:
        destripe_counts(adata, counts_key=counts_key, adjusted_counts_key=adjusted_counts_key)

def expand_labels(adata, labels_key="labels", expanded_labels_key="labels_expanded", algorithm="max_bin_distance", max_bin_distance=2, volume_ratio=4, k=4, subset_pca=True):
    """Expand segmentation labels to nearby bins."""
    adata.obs[expanded_labels_key] = adata.obs[labels_key].values.copy()
    coords = adata.obs[["array_row", "array_col"]].values
    labels = adata.obs[labels_key].values
    
    object_mask = (labels != 0)
    full_reference_inds = np.arange(adata.shape[0])[object_mask]
    full_query_inds = np.arange(adata.shape[0])[~object_mask]
    
    ckd = scipy.spatial.cKDTree(coords[object_mask, :])
    dists, hits = ckd.query(x=coords[~object_mask, :], k=k, workers=-1)
    hits = full_reference_inds[hits]
    calls = labels[hits]
    
    label_values, label_counts = np.unique(labels[object_mask], return_counts=True)
    
    if algorithm == "volume_ratio":
        radii = np.sqrt(label_counts / np.pi)
        label_distances = np.ceil((volume_ratio**(1/3) - 1) * radii)
        label_distance_array = np.zeros(np.max(label_values) + 1)
        label_distance_array[label_values] = label_distances
    elif algorithm == "max_bin_distance":
        label_distance_array = np.full(np.max(label_values) + 1, max_bin_distance)
    else:
        raise ValueError("`algorithm` must be 'max_bin_distance' or 'volume_ratio'")
        
    max_call_distance = label_distance_array[calls]
    dists[dists > max_call_distance] = np.inf
    
    min_per_bin = np.min(dists, axis=1, keepdims=True)
    is_hit = (dists == min_per_bin) & np.isfinite(min_per_bin)
    
    # Case 1: Solitary hit
    clear_mask = (np.sum(is_hit, axis=1) == 1)
    clear_query_inds = full_query_inds[clear_mask]
    clear_query_labels = calls[clear_mask, np.argmin(dists[clear_mask, :], axis=1)]
    adata.obs.loc[adata.obs_names[clear_query_inds], expanded_labels_key] = clear_query_labels
    
    # Case 2: Ambiguous hits, break ties with PCA
    ambiguous_mask = (np.sum(is_hit, axis=1) > 1)
    if np.sum(ambiguous_mask) > 0:
        ambiguous_query_inds = full_query_inds[ambiguous_mask]
        
        if subset_pca:
            involved_bins = np.unique(np.concatenate([hits[ambiguous_mask, :].flatten(), ambiguous_query_inds]))
            pca_subset = sc.pp.pca(np.log1p(adata.X[involved_bins, :]))
            pca = np.zeros((adata.shape[0], pca_subset.shape[1]))
            pca[involved_bins, :] = pca_subset
        else:
            pca = sc.pp.pca(np.log1p(adata.X))
            
        eucl_input = pca[hits[ambiguous_mask, :]] - pca[ambiguous_query_inds, :, None]
        eucl_dists = np.linalg.norm(eucl_input, axis=2)
        eucl_dists[~is_hit[ambiguous_mask, :]] = np.inf
        
        ambiguous_query_labels = calls[ambiguous_mask, np.argmin(eucl_dists, axis=1)]
        adata.obs.loc[adata.obs_names[ambiguous_query_inds], expanded_labels_key] = ambiguous_query_labels

def salvage_secondary_labels(adata, primary_label="labels_he_expanded", secondary_label="labels_gex", labels_key="labels_joint"):
    """Fill unassigned bins from a primary label set with a secondary set."""
    adata.obs[labels_key] = adata.obs[primary_label].copy()
    
    unassigned_mask = adata.obs[primary_label] == 0
    secondary_in_unassigned = adata.obs.loc[unassigned_mask, secondary_label]
    primary_labels_present = set(adata.obs[primary_label].unique())
    
    # Find secondary labels that don't overlap with any primary labels
    secondary_to_take = set(secondary_in_unassigned.unique()) - primary_labels_present
    
    offset = np.max(adata.obs[primary_label])
    mask = unassigned_mask & adata.obs[secondary_label].isin(secondary_to_take)
    adata.obs.loc[mask, labels_key] = adata.obs.loc[mask, secondary_label] + offset
    
    adata.obs[labels_key + "_source"] = "none"
    adata.obs.loc[adata.obs[primary_label] > 0, labels_key + "_source"] = "primary"
    adata.obs.loc[mask, labels_key + "_source"] = "secondary"
    
    if "bin2cell" not in adata.uns:
        adata.uns["bin2cell"] = {}
    adata.uns["bin2cell"].setdefault("secondary_label_offset", {})[labels_key] = offset
    
    print(f"Salvaged {len(secondary_to_take)} secondary labels")

def bin_to_cell(adata, labels_key="labels_expanded", spatial_keys=["spatial"], diameter_scale_factor=None):
    """Collapse bins into single cells based on labels."""
    adata_filt = adata[adata.obs[labels_key] != 0].copy()
    
    cell_to_bin = pd.get_dummies(adata_filt.obs[labels_key], sparse=True)
    cell_names = [str(i) for i in cell_to_bin.columns]
    cell_to_bin_sparse = cell_to_bin.sparse.to_coo().tocsr().T
    
    X = cell_to_bin_sparse.dot(adata_filt.X).tocsr()
    cdata = ad.AnnData(X, var=adata_filt.var)
    cdata.obs_names = cell_names
    cdata.obs['object_id'] = [int(i) for i in cell_names]
    cdata.uns['spatial'] = deepcopy(adata.uns['spatial'])
    
    bin_count = np.asarray(cell_to_bin_sparse.sum(axis=1)).flatten()
    row_means = scipy.sparse.diags(1 / bin_count)
    cell_bin_mean_op = row_means.dot(cell_to_bin_sparse)
    
    cdata.obs['bin_count'] = bin_count
    cdata.obs["array_row"] = cell_bin_mean_op.dot(adata_filt.obs["array_row"].values)
    cdata.obs["array_col"] = cell_bin_mean_op.dot(adata_filt.obs["array_col"].values)
    
    if isinstance(spatial_keys, str):
        spatial_keys = [spatial_keys]
    for key in spatial_keys:
        cdata.obsm[key] = cell_bin_mean_op.dot(adata_filt.obsm[key])
        
    if diameter_scale_factor is None:
        diameter_scale_factor = np.sqrt(np.mean(bin_count))
        
    library = list(cdata.uns['spatial'].keys())[0]
    cdata.uns['spatial'][library]['scalefactors']['spot_diameter_fullres'] *= diameter_scale_factor
    
    source_col = labels_key + "_source"
    if source_col in adata_filt.obs.columns:
        mapping = adata_filt.obs[[labels_key, source_col]].drop_duplicates().astype(str).set_index(labels_key).to_dict()[source_col]
        cdata.obs[source_col] = cdata.obs_names.map(mapping)
        
    return cdata