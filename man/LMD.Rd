% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_lmds.R
\name{LMD}
\alias{LMD}
\title{Calculate Localized Marker Detector Score (LMDS) for each gene}
\usage{
LMD(
  expression,
  feature_space,
  knn = 5,
  kernel = FALSE,
  max_time = 2^20,
  adjust_bridge = TRUE,
  self_loop = 1,
  score_correction = FALSE,
  largeData = TRUE,
  highres = FALSE,
  min_cell = 5
)
}
\arguments{
\item{expression}{matrix; the gene by cell expression matrix.}

\item{feature_space}{matrix; the cell by coordinate matrix (e.g., 20 principal components).}

\item{knn}{integer; the number of nearest neighbors for constructing the graph. Default is 5.}

\item{kernel}{logical; if TRUE, uses a Gaussian kernel. Otherwise, uses a kNN binarized graph. Default is FALSE.}

\item{max_time}{integer; the maximum diffusion time. The actual maximum diffusion time may be shorter if all genes converge beforehand. Default is 2^20.}

\item{adjust_bridge}{logical; if TRUE, connects disconnected components of the graph using Minimum Spanning Trees. Default is TRUE.}

\item{self_loop}{integer; the weight for self-connections. Default is 1.}

\item{score_correction}{logical; if TRUE, adjusts the LMD profile by delta correction. Default is FALSE.}

\item{largeData}{logical; if TRUE, uses functions optimized for large matrix multiplication. Default is TRUE.}

\item{highres}{logical; if TRUE, uses fine time scales (0, 1, 2, 3, ...). If FALSE, uses dyadic time scales (0, 2, 4, ...). Default is FALSE.}

\item{min_cell}{integer; removes genes expressing in fewer than this number of cells. Default is 5.}
}
\value{
A list containing:
\item{score_profile}{data frame; the computed score profile for different diffusion times.}
\item{cumulative_score}{numeric vector; the LMD score for each gene.}
}
\description{
Computes the LMD score of each gene given an expression matrix and cell space.
}
