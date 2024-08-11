# ======================================
# Obtain Modules Activity Score ========
# ======================================
# Extrinsic Function: AddModuleActivityScore

Obtain_cell_partition <- function(expr_dat, gene_partition, 
                                  cell_kNN_graph = NULL, major_vote = 5){
  # Fit GMM
  if(sum(!names(gene_partition) %in% rownames(expr_dat))){
    stop("gene name doesn't match data")}
  
  meta_g = sapply(levels(gene_partition),function(topic){
    if(sum(gene_partition == topic, na.rm = TRUE) > 1){
      return(colMeans(expr_dat[names(gene_partition)[gene_partition == topic],,drop = FALSE]))
    }else{
      return(NULL)
    }
  })
  rownames(meta_g) = colnames(expr_dat)
  
  # tic()
  cell_block = apply(meta_g,2,function(vec){GMM_partition(vec)})
  # toc()
  # cell_block = as.data.frame(meta_g) %>% reframe(across(everything(), GMM_partition))
  cell_block = as.matrix(cell_block) + 0
  
  # Local smooth
  if(!is.null(cell_kNN_graph)){
    cell_block = cell_kNN_graph %*% cell_block
    cell_block = (cell_block >= major_vote) + 0
  }
  
  colnames(cell_block) = levels(gene_partition)
  rownames(cell_block) = colnames(expr_dat)
  return(cell_block)
}
GMM_partition <- function(vector){
  # cat("bottleneck\n")
  i = 1
  opt = ClusterR::Optimal_Clusters_GMM(as.matrix(vector),max_clusters = 10, criterion = "BIC", km_iter = 10, em_iter = 5, seed = i, plot_data = FALSE)
  opt_num = which.min(opt)
  gmm = ClusterR::GMM(as.matrix(vector),gaussian_comps = opt_num, km_iter = 10, em_iter = 5, seed = i)
  labels = stats::predict(gmm,as.matrix(vector))
  
  
  if(length(unique(labels))>2){
    centroid = setNames(gmm$centroids[,1],1:nrow(gmm$centroids))
    
    # Merge each cluster into two groups using hclust
    clusters = cutree(hclust(dist(centroid), method = "average"), k = 2)
    
    labels = labels %in% names(clusters)[clusters == max(clusters)]
  }else{
    labels = labels == which.max(gmm$centroids)
  }
  if(length(table(labels)) > 1){
    if(mean(vector[labels]) < mean(vector[!labels])){
      labels = !labels
    }
  }else{
    if(sum(!labels) == 0){labels = !labels}
  }
  return(labels)
}
GMM_subsampling <- function(seed, gene_partition, expr_dat, cell_kNN_graph, major_vote){
  gene_partition = as.factor(gene_partition)
  set.seed(seed)
  sub_gene_partition = data.frame(gene_partition,
                                  gene_name = names(gene_partition)) %>% 
    group_by(gene_partition) %>% 
    do({
      frac <- case_when(
        length(.$gene_partition) == 1 ~ 1,
        TRUE ~ 0.5
      )
      sample_frac(., frac, replace = FALSE)
    })
  sub_gene_partition = setNames(sub_gene_partition$gene_partition,sub_gene_partition$gene_name)
  return(Obtain_cell_partition(expr_dat, gene_partition = sub_gene_partition, cell_kNN_graph, major_vote))
}

#' Add Module Activity Score
#'
#' Adds module activity scores to a Seurat object based on gene partitions and optionally smooths scores locally.
#'
#' @param srat Seurat object; the single-cell RNA-seq data.
#' @param gene_partition factor; the gene modules.
#' @param assay character; the assay used for calculating module score. Default is "RNA".
#' @param do_local_smooth logical; if TRUE, smooths the module score for each cell by its neighborhoods. Default is FALSE.
#' @param knn integer; the neighborhood size for local smoothing. Default is 10.
#' @param major_vote integer; the majority vote number for local smoothing. Default is 6.
#' @param nloop integer; the number of sampling iterations. Default is 100.
#' @param module_names character vector; the prefix name for each gene module. If NULL, default names are used.
#'
#' @return A Seurat object with module activity scores added to the metadata.
#'
#' @details This function computes module activity scores for gene partitions in a Seurat object. It allows for optional local smoothing of the scores based on neighborhood information.
#'
#' @examples
#' \dontrun{
#'   srat <- AddModuleActivityScore(srat, gene_partition, assay = "RNA", do_local_smooth = TRUE, knn = 10, major_vote = 6, nloop = 100, module_names = NULL)
#' }
#'
#' @importFrom Seurat AddMetaData
#' @export
AddModuleActivityScore <- function(srat, gene_partition, assay = "RNA", 
                                   do_local_smooth = FALSE, knn = 10, major_vote = 6, 
                                   nloop = 100, module_names = NULL){
  # adjust the format of gene_partition
  gene_partition = setNames(as.character(gene_partition),names(gene_partition))
  gene_partition = gene_partition[names(gene_partition) %in% rownames(srat)]
  gene_partition = gene_partition[!is.na(gene_partition)]
  gene_partition = as.factor(gene_partition)
  dat = as.matrix(srat[[assay]]@data[names(gene_partition),])
  if(do_local_smooth){
    if ("pca" %in% names(srat@reductions)){
      ndims = FindPC(srat)
      feature_space <- Embeddings(srat, reduction = "pca")[,1:ndims]
      cat(ncol(feature_space),"PC used for building graph for majority vote\n")
    }else if ("lsi" %in% names(srat@reductions)){
      ndims = 50
      feature_space <- Embeddings(srat, reduction = "lsi")[,2:ndims]
      cat(ncol(feature_space),"LSI used for building graph for majority vote\n")
    }
    else{
      stop("Please RunPCA")
    }
    A = ConstructKnnGraph(knn = knn, feature_space = feature_space)$'adj_matrix'
    cat(sprintf("Do Major Vote: at least %d out of %d neighbors(self-include) expressed\n", major_vote, knn + 1))
  }else{
    A = NULL
  }
  cat("Start Loop\n")
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = nloop,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  cell_block_loop = NULL
  for(loop in 1:nloop){
    cell_block_loop = c(cell_block_loop,list(GMM_subsampling(seed = loop, 
                                                             gene_partition, 
                                                             expr_dat = dat, 
                                                             cell_kNN_graph = A, major_vote)))
    pb$tick()
  }
  cell_block_prop = Reduce(`+`, cell_block_loop) / length(cell_block_loop)
  cell_block = cell_block_prop
  if(is.null(module_names) | length(module_names) != ncol(cell_block)){
    colnames(cell_block) = paste0("Module",colnames(cell_block))
  }else{
    colnames(cell_block) = module_names
  }
  srat <- Seurat::AddMetaData(srat, cell_block, col.name = colnames(cell_block))
  return(srat)
}



