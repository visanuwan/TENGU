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

## Citation

Please cite the following article if you use TENGU in your research
> Visanu Wanchai, Nancy C. Bustamante-Gomez, Alongkorn Kurilung, Karen E. Beenken, Sergio Cortes, Mark S. Smeltzer, Yuet-Kin Leung, Jinhu Xiong, Maria Almeida, Charles A. O'Brien, Intawat Nookaew, A transcriptomic-driven segmentation and cell simulation framework for high-resolution spatial transcriptomics and cell-cell communication, *[Manuscript submitted for publication]* 2026.<br>

---

This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg
