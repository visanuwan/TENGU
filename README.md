# TENGU 👺

TENGU (Transcript-signal ENrichment and Grouping Unit) is a Python package for processing and analyzing 10x Genomics Visium HD spatial transcriptomics data. It provides tools for aggregating 2µm bins into expression profiles of segmented units.

## Features

-   Aggregation of transcription peaks to approximate cell boundaries with a cell simulation option.
-   Deconvolution of bin-level data to `AnnData` or `zarr` objects of segmented units.
-   Secondary analyses such as cell type annotation and cell-cell communication.

## Installation

You can install the package directly from this GitHub repository:

```bash
pip install git+https://github.com/visanuwan/TENGU.git
```

## Testing (Visium HD Mouse Brain)
For this demonstration, we use a publicly available `10X Visium HD Mouse Brain (FFPE)` dataset analyzed using Space Ranger 4.0.1. The input, output, and supplemental files for the Visium HD Gene Expression dataset can be downloaded from [Visium HD Spatial Gene Expression Library, Mouse Brain (FFPE)](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he-v4). The Slide File (H1-8CPZTMY.vlf) required by Space Ranger in this demonstration can be downloaded from [this link](https://www.10xgenomics.com/support/software/space-ranger/2.1/analysis/inputs/image-slide-layout-download). The predictive model used for cell type annotation can be downloaded from [CellTypist's Mouse_Whole_Brain](https://celltypist.cog.sanger.ac.uk/models/Adult_MouseBrain_Yao/v1/Mouse_Whole_Brain.pkl). The cell-cell communication database file can be downloaded from [CellNEST database](https://github.com/schwartzlab-methods/CellNEST/blob/main/database/CellNEST_database.csv).
```
## Go to an example folder
cd example

## Download outputs of 10X Visium HD Mouse Brain and put them into a counts folder
## Download inputs of 10X Visium HD Mouse Brain and put them into a data folder
## Download CellTypist's Mouse_Whole_Brain and put it into a model folder
## Download CellNEST's cell-cell communication database file and put it into a database folder

## (1) Run tengu segmentation
tengu segmentation --counts counts --image data/Visium_HD_Mouse_Brain_tissue_image.tif --out_h5ad TENGU_segmented.h5ad --out_cellmark TENGU_cellmark.csv

## (2) Run Space Ranger with segmented units identified by TENGU
spaceranger count --id=TENGU_Mouse_Brain --fastqs=data/Visium_HD_Mouse_Brain_fastqs --transcriptome=/path/to/refdata-gex-GRCm39-2024-A --probe-set=/path/to/probe_sets/Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv --slide=H1-8CPZTMY --jobmode=local --area=A1 --image=data/Visium_HD_Mouse_Brain_tissue_image.tif --create-bam=true --cytaimage=data/Visium_HD_Mouse_Brain_cytaimage.tif --umi-registration=false --filter-probes=false --custom-segmentation-file=TENGU_cellmark.csv --nucleus-expansion-distance-micron=0 --output-dir=TENGU_counts

## (3) Run tengu annotation
tengu annotation --counts TENGU_counts --model model/Mouse_Whole_Brain.pkl --out predicted_celltype.csv

## (4) Run tengu communication
tengu communication --counts TENGU_counts --database database/CellNEST_database.csv --out ccc.csv

## (5) Run tengu simulation.
## The cellmark file from this step can be used with Space Ranger to generate segmented results, which can be used for TENGU's secondary analyses (annotation and communication).
tengu simulation --counts TENGU_counts --out_zarr TENGU_simulation.zarr --out_cellmark TENGU_simulation_cellmark.csv
```
## Interface
 - `tengu segmentation` - aggregate expression profiles of 2µm bins into segmented units
 - `tengu annotation` - annotate cell types based on a supervised machine-learning classifier
 - `tengu communication` - derive a basic spatially aware cell-cell communication
 - `tengu simulation` - denoise and resolve transcript diffusion with a probabilistic cell simulation framework

## Citation

Please cite the following article if you use TENGU in your research
> Visanu Wanchai, Nancy C. Bustamante-Gomez, Alongkorn Kurilung, Karen E. Beenken, Sergio Cortes, Mark S. Smeltzer, Yuet-Kin Leung, Jinhu Xiong, Maria Almeida, Charles A. O'Brien, Intawat Nookaew, A transcriptomic-driven segmentation and cell simulation framework for high-resolution spatial transcriptomics and cell-cell communication. bioRxiv 2026.04.24.720489; doi: https://doi.org/10.64898/2026.04.24.720489<br>

---

TENGU is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
