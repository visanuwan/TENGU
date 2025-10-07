import numpy as np
import scipy.sparse
from .utils import load_image, normalize

def stardist(image_path, labels_npz_path, stardist_model="2D_versatile_he", block_size=4096, min_overlap=128, context=128, **kwargs):
    """Segment an image with StarDist."""
    from stardist.models import StarDist2D
    
    is_fluo = stardist_model == "2D_versatile_fluo"
    img = load_image(image_path, gray=is_fluo, dtype=np.float16)
    img = normalize(img)
    
    model = StarDist2D.from_pretrained(stardist_model)
    model_axes = "YX" if is_fluo else "YXC"
    
    labels, _ = model.predict_instances_big(
        img, axes=model_axes, block_size=block_size, 
        min_overlap=min_overlap, context=context, **kwargs
    )
    
    labels_sparse = scipy.sparse.csr_matrix(labels)
    scipy.sparse.save_npz(labels_npz_path, labels_sparse)
    print(f"Found {len(np.unique(labels_sparse.data))} objects")