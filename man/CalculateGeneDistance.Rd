% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/obtain_gene_module.R
\name{CalculateGeneDistance}
\alias{CalculateGeneDistance}
\title{Calculate Gene Pairwise Distance}
\usage{
CalculateGeneDistance(dat, method)
}
\arguments{
\item{dat}{matrix; the expression data with genes as rows and cells as columns.}

\item{method}{character; the metric for measuring gene-gene distance. Possible values are "pearson", "euclidean", "KL", "jaccard", and "spearman".
\itemize{
\item "pearson": Calculates Pearson correlation distance (1 - correlation).
\item "euclidean": Calculates Euclidean distance.
\item "KL": Calculates Kullback-Leibler divergence.
\item "jaccard": Calculates Jaccard distance.
\item "spearman": Calculates Spearman correlation distance (1 - correlation).
}}
}
\value{
A distance object representing the pairwise distances between genes.
}
\description{
Calculates the pairwise distance between genes in an expression matrix using various metrics.
}
\examples{
expression_matrix <- matrix(runif(100), nrow=10)
rownames(expression_matrix) <- paste0("Gene", 1:10)
colnames(expression_matrix) <- paste0("Sample", 1:10)
dist <- CalculateGeneDistance(expression_matrix, method="pearson")

}
