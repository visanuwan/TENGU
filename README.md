# TENGU 👺

TENGU (Transcript-signal ENrichment and Grouping Unit) is a Python package for processing and analyzing 10x Genomics Visium HD spatial transcriptomics data. It provides tools for aggregating 2µm bins into cell-level expression profiles.

## Features

-   Aggregation of transcription peaks to approximate cell boundaries.
-   Image segmentation using pre-trained StarDist models.
-   Deconvolution of bin-level data to cell-level `AnnData` objects.

## Installation

You can install the package directly from this GitHub repository:

```bash
pip install git+https://github.com/visanuwan/TENGU.git
```

## Citation

Please cite the following article if you use TENGU in your research
> Visanu Wanchai, Cecile Bustamante-Gomez, Alongkorn Kurilung, Maria Almeida, Charles O’Brien, Intawat Nookaew, TENGU: A Computational Pipeline for High-Fidelity Cellular Deconvolution of Spatial Transcriptomics in Mesenchymal Cells, *[Manuscript in preparation]* 2025.<br>
