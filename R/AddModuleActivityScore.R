#' Add module activity scores
#'
#' @param srat seurat object
#' @param gene_partition vector; gene modules
#' @param assay string; asssay used for calculated modulescore
#' @param do_local_smooth Boolean; TRUE for smooth the module score for each cell by its neighborhoods
#' @param knn integer; neighborhoods size for local smoothing
#' @param major_vote integer; majority vote number for local smoothing
#' @param nloop integer; number of sampling
#' @param module_names string; prefix name for each gene module
#' 
#' 
#'
#' @return
#' A Seurat object
#'
#' @export
#'
#' @examples
#'
#'



AddModuleActivityScore <- function(srat, gene_partition, assay = "RNA", do_local_smooth = FALSE, knn = 10, major_vote = 6, nloop = 100, module_names = NULL){
  # srat: seurat object
  # gene_partition: gene modules
  # assay: asssay used for calculated modulescore
  # do_local_smooth: TRUE for smooth the module score for each cell by its neighborhoods
  # knn: neighborhoods size for local smoothing
  # major_vote: majority vote number for local smoothing
  # nloop: number of sampling
  # module_names: prefix name for each gene module
  
  
  # adjust the format of gene_partition
  gene_partition = setNames(as.character(gene_partition),names(gene_partition))
  gene_partition = gene_partition[names(gene_partition) %in% rownames(srat)]
  gene_partition = gene_partition[!is.na(gene_partition)]
  gene_partition = as.factor(gene_partition)
  dat = srat[[assay]]@data[names(gene_partition),]
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
    A = Symmetric_KNN_graph(knn = knn, feature_space = feature_space)$'adj_matrix'
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
  srat <- AddMetaData(srat, cell_block, col.name = colnames(cell_block))
  return(srat)
}
