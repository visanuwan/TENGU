import numpy as np
import os
import scipy.stats
import skimage
import cv2
from PIL import Image

# setting needed so PIL can load the large TIFFs
Image.MAX_IMAGE_PIXELS = None

# setting needed so cv2 can load the large TIFFs
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(2**40)

def load_image(image_path, gray=False, dtype=np.uint8):
    img = cv2.imread(image_path)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    if gray:
        img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    return img.astype(dtype, copy=False)

def normalize(img):
    eps = 1e-20
    mi = np.percentile(img, 3)
    ma = np.percentile(img, 99.8)
    return (img - mi) / (ma - mi + eps)

def check_array_coordinates(adata, row_max=3349, col_max=3349):
    if "bin2cell" not in adata.uns:
        adata.uns["bin2cell"] = {}
    adata.uns["bin2cell"]["array_check"] = {}
    for axis in ["row", "col"]:
        adata.uns["bin2cell"]["array_check"][axis] = {}
        if axis == "row":
            adata.uns["bin2cell"]["array_check"][axis]["max"] = row_max
        elif axis == "col":
            adata.uns["bin2cell"]["array_check"][axis]["max"] = col_max
        
        single_axis, spatial_axis = ("row", 0) if axis == "col" else ("col", 1)
        
        val = adata.obs['array_' + single_axis].value_counts().index[0]
        mask = (adata.obs['array_' + single_axis] == val)
        array_vals = adata.obs.loc[mask, 'array_' + axis].values
        spatial_vals = adata.obsm['spatial'][mask, spatial_axis]
        
        adata.uns["bin2cell"]["array_check"][axis]["flipped"] = scipy.stats.pearsonr(array_vals, spatial_vals)[0] < 0

def actual_vs_inferred_image_shape(adata, img, ratio_threshold=0.99):
    library = list(adata.uns['spatial'].keys())[0]
    inferred_dim = np.array(adata.uns['spatial'][library]['images']['hires'].shape)[:2] / adata.uns['spatial'][library]['scalefactors']['tissue_hires_scalef']
    actual_dim = np.array(img.shape)[:2]
    if np.min(np.hstack((actual_dim / inferred_dim, inferred_dim / actual_dim))) < ratio_threshold:
        raise ValueError(f"Morphology image dimension mismatch. Inferred: {inferred_dim}, Actual: {actual_dim}. Ensure --image in Space Ranger matches source_image_path.")

def mpp_to_scalef(adata, mpp):
    library = list(adata.uns['spatial'].keys())[0]
    mpp_source = adata.uns['spatial'][library]['scalefactors']['microns_per_pixel']
    return mpp_source / mpp

def get_mpp_coords(adata, basis="spatial", spatial_key="spatial", mpp=None):
    if basis == "array" and mpp is None:
        raise ValueError("Must specify mpp for array coordinates.")
    
    if basis == "spatial":
        scalef = mpp_to_scalef(adata, mpp=mpp) if mpp is not None else 1
        coords = (adata.obsm[spatial_key] * scalef).astype(int)[:, ::-1]
    elif basis == "array":
        scalef = 2 / mpp
        coords = np.round(adata.obs[['array_row', 'array_col']].values * scalef).astype(int)
        if adata.uns["bin2cell"]["array_check"]["row"]["flipped"]:
            coords[:, 0] = np.round(adata.uns["bin2cell"]["array_check"]["row"]["max"] * scalef).astype(int) - coords[:, 0]
        if adata.uns["bin2cell"]["array_check"]["col"]["flipped"]:
            coords[:, 1] = np.round(adata.uns["bin2cell"]["array_check"]["col"]["max"] * scalef).astype(int) - coords[:, 1]
    return coords

def get_crop(adata, basis="spatial", spatial_key="spatial", mpp=None, buffer=0):
    coords = get_mpp_coords(adata, basis=basis, spatial_key=spatial_key, mpp=mpp)
    return (
        np.max([np.min(coords[:, 1]) - buffer, 0]),
        np.max([np.min(coords[:, 0]) - buffer, 0]),
        np.max(coords[:, 1]) + buffer,
        np.max(coords[:, 0]) + buffer
    )

def _create_scaled_image(adata, img_path, is_he, channel, mpp, crop, buffer, spatial_cropped_key, store, img_key, save_path):
    import tifffile as tf
    library = list(adata.uns['spatial'].keys())[0]
    
    if is_he:
        img = load_image(img_path)
    else: # is_if
        img = tf.imread(img_path, key=channel).astype(np.float32)
        img = normalize(img)
        img[img < 0] = 0
        img[img > 1] = 1

    actual_vs_inferred_image_shape(adata, img)

    if crop:
        crop_coords = get_crop(adata, basis="spatial", mpp=None, buffer=buffer)
        img = img[crop_coords[1]:crop_coords[3], crop_coords[0]:crop_coords[2]]
        if spatial_cropped_key is None:
            spatial_cropped_key = f"spatial_cropped_{buffer}_buffer"
        adata.obsm[spatial_cropped_key] = adata.obsm["spatial"].copy()
        adata.obsm[spatial_cropped_key][:, 0] -= crop_coords[0]
        adata.obsm[spatial_cropped_key][:, 1] -= crop_coords[1]
        print(f"Cropped spatial coordinates key: {spatial_cropped_key}")

    scalef = mpp_to_scalef(adata, mpp=mpp)
    dim = (np.array(img.shape[:2]) * scalef).astype(int)[::-1]
    img = cv2.resize(img, dim, interpolation=cv2.INTER_AREA)

    if store:
        if img_key is None:
            img_key = f"{mpp}_mpp"
            if crop:
                img_key += f"_{buffer}_buffer"
        adata.uns['spatial'][library]['images'][img_key] = img
        adata.uns['spatial'][library]['scalefactors'][f'tissue_{img_key}_scalef'] = scalef
        print(f"Image key: {img_key}")

    if save_path is not None:
        if is_he:
            cv2.imwrite(save_path, cv2.cvtColor(img, cv2.COLOR_RGB2BGR))
        else:
            cv2.imwrite(save_path, cv2.cvtColor((255 * img).astype(np.uint8), cv2.COLOR_GRAY2BGR))

def scaled_he_image(adata, mpp=1, crop=True, buffer=150, **kwargs):
    library = list(adata.uns['spatial'].keys())[0]
    img_path = adata.uns['spatial'][library]['metadata']['source_image_path']
    _create_scaled_image(adata, img_path, is_he=True, channel=None, mpp=mpp, crop=crop, buffer=buffer, **kwargs)

def scaled_if_image(adata, channel, mpp=1, crop=True, buffer=150, **kwargs):
    library = list(adata.uns['spatial'].keys())[0]
    img_path = adata.uns['spatial'][library]['metadata']['source_image_path']
    _create_scaled_image(adata, img_path, is_he=False, channel=channel, mpp=mpp, crop=crop, buffer=buffer, **kwargs)

def insert_labels(adata, labels_npz_path, basis="spatial", spatial_key="spatial", mpp=None, labels_key="labels"):
    import scipy.sparse
    labels_sparse = scipy.sparse.load_npz(labels_npz_path)
    
    if "bin2cell" not in adata.uns:
        adata.uns["bin2cell"] = {}
    if "labels_npz_paths" not in adata.uns["bin2cell"]:
        adata.uns["bin2cell"]["labels_npz_paths"] = {}
        
    adata.uns["bin2cell"]["labels_npz_paths"][labels_key] = os.path.abspath(labels_npz_path)
    
    coords = get_mpp_coords(adata, basis=basis, spatial_key=spatial_key, mpp=mpp)
    
    adata.obs[labels_key] = 0
    mask = (
        (coords[:, 0] >= 0) & (coords[:, 0] < labels_sparse.shape[0]) &
        (coords[:, 1] >= 0) & (coords[:, 1] < labels_sparse.shape[1])
    )
    adata.obs.loc[mask, labels_key] = np.asarray(labels_sparse[coords[mask, 0], coords[mask, 1]]).flatten()

def grid_image(adata, val, log1p=False, mpp=2, sigma=None, save_path=None):
    import scipy.sparse
    if val in adata.obs.columns:
        vals = adata.obs[val].values.copy()
    elif val in adata.var_names:
        vals = adata[:, val].X
        if scipy.sparse.issparse(vals):
            vals = vals.toarray()
        vals = np.asarray(vals).flatten()
    else:
        raise ValueError(f'"{val}" not found in .obs or .var_names')

    vals = (255 * (vals - np.min(vals)) / (np.max(vals) - np.min(vals))).astype(np.uint8)
    if log1p:
        vals = np.log1p(vals)
        vals = (255 * (vals - np.min(vals)) / (np.max(vals) - np.min(vals))).astype(np.uint8)

    if "bin2cell" not in adata.uns or "array_check" not in adata.uns["bin2cell"]:
        check_array_coordinates(adata)

    img = np.zeros((adata.uns["bin2cell"]["array_check"]["row"]["max"] + 1,
                    adata.uns["bin2cell"]["array_check"]["col"]["max"] + 1), dtype=np.uint8)
    img[adata.obs['array_row'], adata.obs['array_col']] = vals

    if adata.uns["bin2cell"]["array_check"]["row"]["flipped"]:
        img = np.flip(img, axis=0)
    if adata.uns["bin2cell"]["array_check"]["col"]["flipped"]:
        img = np.flip(img, axis=1)

    if mpp != 2:
        dim = np.round(np.array(img.shape) * 2 / mpp).astype(int)[::-1]
        img = cv2.resize(img, dim, interpolation=cv2.INTER_CUBIC)

    if sigma is not None:
        img = skimage.filters.gaussian(img, sigma=sigma)
        img = (255 * (img - np.min(img)) / (np.max(img) - np.min(img))).astype(np.uint8)

    if save_path is not None:
        cv2.imwrite(save_path, img)
    else:
        return img