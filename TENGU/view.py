import numpy as np
import scipy.sparse
import scipy.stats
import skimage.segmentation
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from PIL import Image
from .utils import load_image, normalize

def view_labels(image_path, labels_npz_path, crop=None, stardist_normalize=False, fill=False, border=True, fill_palette=None, fill_label_weight=0.5, border_color=[255, 255, 0]):
    """Render segmentation results on a full or cropped image."""
    labels_sparse = scipy.sparse.load_npz(labels_npz_path)
    dtype = np.float16 if stardist_normalize else np.uint8

    if crop is None:
        img = load_image(image_path, dtype=dtype)
    else:
        pil_img = Image.open(image_path)
        img = np.array(pil_img.crop(crop).convert('RGB'), dtype=dtype)
        labels_sparse = labels_sparse[crop[1]:crop[3], crop[0]:crop[2]]

    if stardist_normalize:
        img = normalize(img)
        img = np.clip(img, 0, 1)
        img = (255 * img).astype(np.uint8)

    labels_coo = labels_sparse.tocoo()

    if fill:
        if fill_palette is None:
            palette = (np.array(sns.color_palette("bright")) * 255).astype(np.uint8)
            fill_palette = np.delete(palette, 6, 0)
        
        img[labels_coo.row, labels_coo.col, :] = \
            (1 - fill_label_weight) * img[labels_coo.row, labels_coo.col, :] + \
            fill_label_weight * fill_palette[labels_coo.data % fill_palette.shape[0], :]

    if border:
        border_sparse = scipy.sparse.coo_matrix(skimage.segmentation.find_boundaries(np.array(labels_sparse.todense())))
        img[border_sparse.row, border_sparse.col, :] = border_color
        
    return img

def _overlay_onto_img(img, labels_sparse, cdata, key, common_objects, fill_label_weight=1, cat_cmap="tab20", cont_cmap="viridis"):
    """Helper function to overlay metadata/expression onto an image."""
    fig, ax = plt.subplots(figsize=(5, 2))
    ax.axis("off")

    if key in cdata.obs.columns and ("category" in cdata.obs[key].dtype.name):
        # Categorical data
        cats = np.zeros(np.max(common_objects) + 1, dtype=np.int32)
        cats[common_objects] = cdata.obs.loc[[str(i) for i in common_objects], key].cat.codes
        cats_unique_original = np.unique(cats[common_objects])
        cats[common_objects] = scipy.stats.rankdata(cats[common_objects], method="dense") - 1
        
        fill_palette = (np.array(sns.color_palette(cat_cmap, n_colors=len(cats_unique_original))) * 255).astype(np.uint8)
        
        img[labels_sparse.row, labels_sparse.col, :] = \
            (1 - fill_label_weight) * img[labels_sparse.row, labels_sparse.col, :] + \
            fill_label_weight * fill_palette[cats[labels_sparse.data], :]
        
        legend_patches = [
            matplotlib.patches.Patch(color=fill_palette[i, :] / 255.0, label=cdata.obs[key].cat.categories[j])
            for i, j in zip(np.unique(cats[common_objects]), cats_unique_original)
        ]
        ax.legend(handles=legend_patches, loc="center", title=key, frameon=False)

    else: # Continuous data (from .obs or .var)
        if key in cdata.obs.columns:
            vals = cdata.obs.loc[[str(i) for i in common_objects], key].values
        elif key in cdata.var_names:
            vals = cdata[[str(i) for i in common_objects]][:, key].X.toarray().flatten()
        else:
            raise ValueError(f"'{key}' not found in cdata.obs or cdata.var_names")
        
        colormap = matplotlib.colormaps.get_cmap(cont_cmap)
        norm = matplotlib.colors.Normalize(vmin=np.min(vals), vmax=np.max(vals))
        fig.colorbar(matplotlib.cm.ScalarMappable(cmap=colormap, norm=norm), ax=ax, orientation="horizontal", label=key)
        
        scaled_vals = (vals - np.min(vals)) / (np.max(vals) - np.min(vals) + 1e-9)
        fill_palette = np.zeros((np.max(common_objects) + 1, 3))
        fill_palette[common_objects, :] = (colormap(scaled_vals)[:, :3] * 255).astype(np.uint8)
        
        img[labels_sparse.row, labels_sparse.col, :] = \
            (1 - fill_label_weight) * img[labels_sparse.row, labels_sparse.col, :] + \
            fill_label_weight * fill_palette[labels_sparse.data, :]

    plt.close(fig)
    return img, fig

def view_cell_labels(image_path, labels_npz_path, cdata, fill_key=None, border_key=None, crop=None, stardist_normalize=False, fill_label_weight=1, thicken_border=True, **kwargs):
    """Color morphology segmentations by cell-level metadata or gene expression."""
    labels_sparse = scipy.sparse.load_npz(labels_npz_path)
    dtype = np.float16 if stardist_normalize else np.uint8

    if crop is None:
        img = load_image(image_path, dtype=dtype)
    else:
        pil_img = Image.open(image_path)
        img = np.array(pil_img.crop(crop).convert('RGB'), dtype=dtype)
        labels_sparse = labels_sparse[crop[1]:crop[3], crop[0]:crop[2]]

    if stardist_normalize:
        img = normalize(img)
        img = np.clip(img, 0, 1)
        img = (255 * img).astype(np.uint8)

    labels_coo = labels_sparse.tocoo()
    common_objects = np.sort(list(set(np.unique(labels_coo.data)).intersection(set([int(i) for i in cdata.obs_names]))))
    
    labels_coo.data[~np.isin(labels_coo.data, common_objects)] = 0
    labels_coo.eliminate_zeros()
    
    legends = {}

    if fill_key is not None:
        img, fig = _overlay_onto_img(img.copy(), labels_coo, cdata, fill_key, common_objects, fill_label_weight, **kwargs)
        legends[fill_key] = fig

    if border_key is not None:
        border = skimage.segmentation.find_boundaries(np.array(labels_sparse.todense()), mode="inner")
        coords = np.nonzero(border)
        if thicken_border:
            border_rows = np.hstack([np.clip(coords[0] + i, 0, border.shape[0] - 1) for i in [-1, 1, 0, 0]])
            border_cols = np.hstack([np.clip(coords[1] + i, 0, border.shape[1] - 1) for i in [0, 0, -1, 1]])
            border[border_rows, border_cols] = True
            coords = np.nonzero(border)
        
        border_labels_coo = scipy.sparse.coo_matrix(
            (np.array(labels_sparse.tocsr()[coords]).flatten(), coords), 
            shape=labels_sparse.shape
        )
        border_labels_coo.eliminate_zeros()
        
        img, fig = _overlay_onto_img(img, border_labels_coo, cdata, border_key, common_objects, fill_label_weight=1, **kwargs)
        legends[border_key] = fig
        
    return img, legends