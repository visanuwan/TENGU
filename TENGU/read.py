from pathlib import Path
import json
import pandas as pd
from matplotlib.image import imread
from scanpy import read_10x_h5
from scanpy import logging as logg
import anndata as ad

def read_visium(
    path: Path | str,
    genome: str | None = None,
    *,
    count_file: str = "filtered_feature_bc_matrix.h5",
    library_id: str | None = None,
    load_images: bool | None = True,
    source_image_path: Path | str | None = None,
    spaceranger_image_path: Path | str | None = None,
) -> ad.AnnData:
    """Reads 10x Genomics Visium data."""
    path = Path(path)
    if spaceranger_image_path is None:
        spaceranger_image_path = path / "spatial"
    else:
        spaceranger_image_path = Path(spaceranger_image_path)
    
    adata = read_10x_h5(path / count_file, genome=genome)
    adata.uns["spatial"] = {}

    from h5py import File
    with File(path / count_file, mode="r") as f:
        attrs = dict(f.attrs)
    if library_id is None:
        library_id = str(attrs.pop("library_ids")[0], "utf-8")

    adata.uns["spatial"][library_id] = {}

    if load_images:
        positions_path = path / "spatial"
        pos_files = ["tissue_positions.csv", "tissue_positions.parquet", "tissue_positions_list.csv"]
        tissue_positions_file = next((positions_path / f for f in pos_files if (positions_path / f).exists()), None)

        if not tissue_positions_file:
            raise FileNotFoundError("Could not find tissue positions file.")

        files = {
            "scalefactors_json_file": path / "spatial/scalefactors_json.json",
            "hires_image": spaceranger_image_path / "tissue_hires_image.png",
            "lowres_image": spaceranger_image_path / "tissue_lowres_image.png",
        }

        for f_path in files.values():
            if not f_path.exists():
                if "image" in str(f_path):
                    logg.warning(f"Could not find image file: '{f_path}'")
                else:
                    raise OSError(f"Could not find required file: '{f_path}'")

        adata.uns["spatial"][library_id]["images"] = {}
        for res in ["hires", "lowres"]:
            try:
                adata.uns["spatial"][library_id]["images"][res] = imread(str(files[f"{res}_image"]))
            except FileNotFoundError:
                 logg.warning(f"Could not find '{res}_image'.")
            except Exception as e:
                raise OSError(f"Error reading '{res}_image': {e}")


        adata.uns["spatial"][library_id]["scalefactors"] = json.loads(files["scalefactors_json_file"].read_bytes())
        adata.uns["spatial"][library_id]["metadata"] = {
            k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
            for k in ("chemistry_description", "software_version") if k in attrs
        }

        if tissue_positions_file.name.endswith(".csv"):
            positions = pd.read_csv(tissue_positions_file, header=0 if tissue_positions_file.name == "tissue_positions.csv" else None, index_col=0)
        elif tissue_positions_file.name.endswith(".parquet"):
            positions = pd.read_parquet(tissue_positions_file)
            positions.set_index("barcode", inplace=True)
        
        positions.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
        adata.obs = adata.obs.join(positions, how="left")
        adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()
        adata.obs.drop(columns=["pxl_row_in_fullres", "pxl_col_in_fullres"], inplace=True)

        if source_image_path is not None:
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(Path(source_image_path).resolve())

    return adata