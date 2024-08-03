#' Construct a symmetric KNN graph
#'
#' @param knn integer; knn for constructing the graph
#' @param feature_space matrix; cell x coordinate (e.g. 20PCs) matrix
#' @param adjust_by_MST boolean; TRUE for connecting disconnected components by Minimum Spanning Trees
#' @param self-loop integer; weight for self connect (default 1)
#'
#'
#' @return
#' A list of objects
#'
#' @export
#'
#' @examples
#'
#'


Symmetric_KNN_graph <- function(knn = 5, feature_space, adjust_by_MST = TRUE, self_loop = 1){
  # knn: knn for constructing the graph
  # feature_space: cell x coordinate (e.g. 20PCs) matrix
  # adjust_by_MST: TRUE for connecting disconnected components by Minimum Spanning Trees
  # self-loop: weight for self connect (default 1)
  
  cat("Constructing KNN graph\n")
  knn_list <- FNN::get.knn(feature_space, k = knn, algorithm = "kd_tree")
  Idx = knn_list$nn.index
  A = Matrix::sparseMatrix(i = rep(1:nrow(Idx),each = knn),
                           j = c(t(Idx)),
                           dims = c(nrow(Idx), nrow(Idx)),
                           x = 1)
  rownames(A) = colnames(A) = rownames(feature_space)
  
  # remove orphans/singletons
  res = findDisconnectedComponents(A)
  filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 10) comp else NULL)
  filtered_components <- Filter(Negate(is.null), filtered_components)
  filtered_node_names <- unlist(lapply(filtered_components, function(comp) names(comp)))
  A = A[rownames(A)%in%filtered_node_names,
        colnames(A)%in%filtered_node_names]
  
  # Connect Components Using MST Principles
  if(adjust_by_MST){
    while(length(filtered_components) > 1){
      cat(length(filtered_components), " Disconnected Components, adjust by MST.\n")
      edgesToAdd <- lapply(1:(length(filtered_components)-1), function(i) {
        lapply((i + 1):length(filtered_components), function(j) {
          findShortestEdge(filtered_components[[i]], filtered_components[[j]], feature_space)
        })
      })
      edgesToAdd = do.call(rbind.data.frame, lapply(unlist(edgesToAdd, recursive = FALSE),function(x) names(x))) %>% distinct()
      for (i in 1:nrow(edgesToAdd)) {
        A[edgesToAdd[i,1],edgesToAdd[i,2]] = 1
      }
      res = findDisconnectedComponents(A)
      filtered_components = res$components
    }
  }
  
  # if(adjust){ # increase knn if containing disconnected components
  #   while(length(filtered_components) >= 2){
  #     knn = knn + 1
  #     knn_list <- FNN::get.knn(feature_space, k = knn, algorithm = "kd_tree")
  #     Idx = knn_list$nn.index
  #     # Adjacency matrix
  #     A = Matrix::sparseMatrix(i = rep(1:nrow(Idx),each = knn),
  #                              j = c(t(Idx)),
  #                              dims = c(nrow(Idx), nrow(Idx)),
  #                              x = 1)
  #     rownames(A) = colnames(A) = rownames(feature_space)
  #     res = findDisconnectedComponents(A)
  #     filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 5) comp else NULL)
  #     filtered_components <- Filter(Negate(is.null), filtered_components)
  #   }
  #   print(paste0("Final kNN: ",knn))
  # }
  
  if (length(filtered_components) >= 2) {
    stop("Disconnected Components, Please increase knn or adjust by MST.")
  }
  
  A = as.matrix(A)
  diag(A) = self_loop
  # Graph affinity matrix by symmetrize the adjacency mtx
  # W = (A + t(A))/2 # weighted
  W = pmin(A + t(A),1) # unweighted
  diag(W) = self_loop
  if(dim(A)[1] < nrow(feature_space)){
    cat("Remove ", nrow(feature_space) - dim(A)[1], " disconnected cells.\n")
  }
  if(ncol(W) > 1e4){
    return(list(graph = as(W,"sparseMatrix"), adj_matrix = as(A,"sparseMatrix"), component = filtered_components))
  }else{
    return(list(graph = W, adj_matrix = A, component = filtered_components))
  }
}
