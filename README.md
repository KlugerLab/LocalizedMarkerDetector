
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Localized Marker Detector

<!-- badges: start -->
<!-- badges: end -->

Localized Marker Detector (LMD) is a computational framework designed
for the identification of gene expression markers localized to specific
cell populations within single-cell RNA sequencing data. The major
workflow of LMD comprises the following three main steps:

-   Step1. Constructing a cell-cell affinity graph
-   Step2. Diffusing the gene expression value across the cell graph
-   Step3. Assigning a score to each gene based on the dynamics of its
    diffusion process
-   Optional Downstream tasks
    -   Identifying gene modules and characterizing functional cell
        groups
    -   Cross-sample comparison

<img src="./man/figures/LMD_workflow.png" width="50%" height="auto" />

## Installation

LMD can be installed in R as follows:

``` r
install.packages("devtools")
devtools::install_github("ruiqi0130/LocalizedMarkerDetector")
```

## Example tutorial

Please check [LMD
tutorial](https://ruiqi0130.github.io/LocalizedMarkerDetector/articles/).

## References

References of LMD functions can be found
[here](https://ruiqi0130.github.io/LocalizedMarkerDetector/reference/index.html).

<!-- Data used in the tutorial can be downloaded from [Figshare](https://figshare.com/articles/dataset/Processed_Seurat_objects_for_GeneTrajectory_inference_Gene_Trajectory_Inference_for_Single-cell_Data_by_Optimal_Transport_Metrics_/25243225). -->
