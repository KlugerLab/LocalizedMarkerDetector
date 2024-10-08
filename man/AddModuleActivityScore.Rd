% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_module_score.R
\name{AddModuleActivityScore}
\alias{AddModuleActivityScore}
\title{Add Module Activity Score}
\usage{
AddModuleActivityScore(
  srat,
  gene_partition,
  assay = "RNA",
  do_local_smooth = FALSE,
  knn = 10,
  major_vote = 6,
  nloop = 100,
  module_names = NULL
)
}
\arguments{
\item{srat}{Seurat object; the single-cell RNA-seq data.}

\item{gene_partition}{factor; the gene modules.}

\item{assay}{character; the assay used for calculating module score. Default is "RNA".}

\item{do_local_smooth}{logical; if TRUE, smooths the module score for each cell by its neighborhoods. Default is FALSE.}

\item{knn}{integer; the neighborhood size for local smoothing. Default is 10.}

\item{major_vote}{integer; the majority vote number for local smoothing. Default is 6.}

\item{nloop}{integer; the number of sampling iterations. Default is 100.}

\item{module_names}{character vector; the prefix name for each gene module. If NULL, default names are used.}
}
\value{
A Seurat object with module activity scores added to the metadata.
}
\description{
Adds module activity scores to a Seurat object based on gene partitions and optionally smooths scores locally.
}
\details{
This function computes module activity scores for gene partitions in a Seurat object. It allows for optional local smoothing of the scores based on neighborhood information.
}
