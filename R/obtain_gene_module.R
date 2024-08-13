# ============================
# Obtain Gene Modules ========
# ============================
# Extrinsic Function: CalculateGeneDistance, ClusterGenes

#' Calculate Gene Pairwise Distance
#'
#' Calculates the pairwise distance between genes in an expression matrix using various metrics.
#'
#' @param dat matrix; the expression data with genes as rows and cells as columns.
#' @param method character; the metric for measuring gene-gene distance. Possible values are "pearson", "euclidean", "KL", "jaccard", and "spearman".
#' \itemize{
#'   \item "pearson": Calculates Pearson correlation distance (1 - correlation).
#'   \item "euclidean": Calculates Euclidean distance.
#'   \item "KL": Calculates Kullback-Leibler divergence.
#'   \item "jaccard": Calculates Jaccard distance.
#'   \item "spearman": Calculates Spearman correlation distance (1 - correlation).
#' }
#' @return A distance object representing the pairwise distances between genes.
#'
#' @examples
#' expression_matrix <- matrix(runif(100), nrow=10)
#' rownames(expression_matrix) <- paste0("Gene", 1:10)
#' colnames(expression_matrix) <- paste0("Sample", 1:10)
#' dist <- CalculateGeneDistance(expression_matrix, method="pearson")
#'
#' @importFrom stats dist
#' @importFrom philentropy distance
#' @export
CalculateGeneDistance <- function(dat, method){
  # dat: expression matrix
  # method: metric for measuring gene-gene distance (option: pearson, euclidean, KL, jaccard, spearman)
  rho_dat = RowwiseNormalize(dat)
  if(method == "pearson"){
    cor = cor(t(rho_dat), use = "pairwise.complete.obs", method = "pearson")
    dist = 1 - cor
  }else if(method == "euclidean"){
    dist = stats::dist(rho_dat, method = "euclidean") 
  }else if(method == "KL"){
    dist = philentropy::distance(rho_dat, method = "kullback-leibler", unit = "log", epsilon = 1e-9)
  }else if(method == "jaccard"){
    dist = philentropy::distance(rho_dat, method = "jaccard")
  }else if(method == "spearman"){
    cor = cor(t(rho_dat), use = "pairwise.complete.obs", method = "spearman")
    dist = 1 - cor
  }
  
  rownames(dist) = colnames(dist) = rownames(rho_dat)
  dist = as.dist(dist)
  return(dist)
}


#' Cluster Genes
#'
#' This function partitions genes based on their pairwise distances using various clustering methods.
#'
#' @param dist dist; the gene-gene distance matrix.
#' @param clustering_method character; the method for clustering the genes. Options are "average", "ward.D", "ward.D2", "single", "complete", "mcquitty", "median", "centroid", "dbscan", or "hdbscan". Default is "average".
#' @param min_gene integer; the minimum number of genes each group should contain. Default is 10.
#' @param deepSplit integer; parameters for \code{cutreeDynamic}. Default is 2.
#' @param return_tree logical; if TRUE, returns a gene partition tree; otherwise returns only the gene partition. Default is FALSE.
#' @param filtered logical; if TRUE, filters out some genes for each partition based on the SML method (Parisi et al., 2014). Default is TRUE.
#' @param accu numeric; the threshold for filtering out genes. Default is 0.75.
#'
#' @return If \code{return_tree} is TRUE, returns a list containing the gene partition and the gene partition tree. Otherwise, returns the gene partition.
#'
#' @details This function performs hierarchical clustering using the specified method and partitions the genes. It also provides options for other clustering methods like dbscan and hdbscan, and to filter out noisy genes.
#'
#' @examples
#' gene_dist <- dist(matrix(runif(100), nrow=10))
#' gene_partition <- ClusterGenes(gene_dist, clustering_method="average")
#' 
#'
#' @importFrom stats hclust
#' @importFrom dynamicTreeCut cutreeDynamic
#' @importFrom dbscan dbscan hdbscan
#' @export
#'
#' @references
#' Parisi, F., Strino, F., Nadler, B., & Kluger, Y. (2014). Ranking and combining multiple predictors without labeled data. Proceedings of the National Academy of Sciences, 111(4), 1253-1258. doi:10.1073/pnas.1219097111
ClusterGenes <- function(dist, clustering_method = "average", 
                                  min_gene = 10, deepSplit = 2, 
                                  return_tree = FALSE, 
                                  filtered = TRUE, accu = 0.75){
  # check if perform hclust
  if(clustering_method %in% c("ward.D","ward.D2","single",
                              "complete","average","mcquitty",
                              "median","centroid")){
    gene_hree = hclust(dist, method = clustering_method)
    gene_partition = dynamicTreeCut::cutreeDynamic(
      dendro = gene_hree, 
      distM = as.matrix(dist), deepSplit = deepSplit,
      pamStage = TRUE, pamRespectsDendro = TRUE,
      minClusterSize = min_gene)
    names(gene_partition) = labels(dist)
    gene_partition = as.factor(gene_partition)
    # Re-name gene modules based on the order of htree
    module_order = order(sapply(levels(gene_partition), function(mod) {
      median(which(gene_partition[gene_hree$order] == mod))
    }))
    levels(gene_partition)[module_order] = seq(nlevels(gene_partition))
    gene_partition = factor(gene_partition,levels = seq(nlevels(gene_partition)))
    if(filtered){
      gene_partition_filter = unlist(lapply(levels(gene_partition),function(i){
        genes = names(gene_partition)[gene_partition == i]
        if(length(genes) < 10){return(NULL)}
        eig_res = RSpectra::eigs_sym((1 - as.matrix(dist))[genes,genes],k=1)
        norm_eig_vec = sqrt(eig_res$values) * abs(eig_res$vectors[,1])
        names(norm_eig_vec) = genes
        norm_eig_vec = sort(norm_eig_vec)
        # p = ggplot(data.frame(idx = 1:length(norm_eig_vec),norm_eigen_vec = norm_eig_vec),aes(x = idx,y = norm_eig_vec)) + geom_point() + geom_hline(yintercept = 0.5) + ylim(0,1) + theme_minimal() + ggtitle(paste0("Module",i))
        genes_left = names(norm_eig_vec)[norm_eig_vec > 2 * accu - 1]
        if(length(genes_left) < 10){return(NULL)}
        else{return(genes_left)}
      }))
      cat("Filtering out outlier genes in each module: ", length(gene_partition_filter)," genes left.\n")
      gene_partition[!names(gene_partition) %in% gene_partition_filter] = NA
      gene_partition = gene_partition[!is.na(gene_partition)]
      gene_partition <- droplevels(gene_partition)
      # gene_partition <- addNA(gene_partition, ifany = TRUE)
    }
    if(return_tree){
      return(list(gene_partition = gene_partition,gene_hree = gene_hree))
    }
    else{return(gene_partition)}
  }
  else if(clustering_method == "dbscan"){
    set.seed(123) # For reproducibility
    db <- dbscan::dbscan(dist, eps = 0.5)
    # dbscan::kNNdistplot(dist, k =  20)
    names(db$cluster) = labels(dist)
    return(db$cluster)
  }
  else if(clustering_method == "hdbscan"){
    db <- dbscan::hdbscan(dist, minPts = 5)
    names(db$cluster) = labels(dist)
    return(db$cluster)
  }
  else{
    stop("Method not found!")
  }
}
