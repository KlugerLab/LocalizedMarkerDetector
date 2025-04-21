# =========================================
# Construct Cell Graph ====================
# =========================================
# Extrinsic Function: ConstructKnnGraph, ConstructGaussianGraph, ConstructDiffusionOperators

findDisconnectedComponents <- function(adj_matrix) {
  # Create a graph from the adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  # Find the connected components
  components <- igraph::components(g)
  # Get the number of disconnected components
  num_components <- components$no
  # Get the corresponding nodes for each component
  component_nodes <- split(igraph::V(g), components$membership)
  list(number_of_components = num_components, components = component_nodes)
}
findShortestEdge <- function(component1, component2, data) {
  # Calculate all pairwise distances between nodes in the two components
  distances <- outer(component1, component2, Vectorize(function(x, y) dist(data[c(x,y), ])))
  # Find the minimum distance and the corresponding nodes
  min_distance_idx <- which(distances == min(distances), arr.ind = TRUE)
  return(list(from = component1[min_distance_idx[1]], to = component2[min_distance_idx[2]], distance = min(distances)))
}

#' Construct a Symmetric KNN graph
#'
#' Constructs a symmetric KNN graph from a given feature space.
#'
#' @param knn integer; the number of nearest neighbors to consider for constructing the graph. Default is 5.
#' @param feature_space matrix; a cell-by-coordinate matrix (e.g., 20 principal components).
#' @param adjust_disconnection boolean; TRUE for fully connecting disconnected components. Default is TRUE.
#' @param self_loop integer; weight for self connections (default is 1).
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{graph}{matrix; symmetric KNN graph W, computed as \code{pmax(1, (A + A^T) / 2)}.}
#'   \item{adj_matrix}{matrix; adjacency matrix A.}
#'   \item{component}{list; disconnected components in the graph, each component is represented as a subgraph.}
#' }
#'
#' @examples
#' feature_space <- matrix(rnorm(200), nrow = 10, ncol = 20)
#' result <- ConstructKnnGraph(knn = 3, feature_space = feature_space)
#' graph <- result$graph
#' adj_matrix <- result$adj_matrix
#' 
#' @export
ConstructKnnGraph <- function(knn = 5, feature_space, 
                              adjust_disconnection = TRUE, self_loop = 1){
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
  if(dim(A)[1] < nrow(feature_space)){
    cat("Remove ", nrow(feature_space) - dim(A)[1], " singleton cells.\n")
  }
  
  # Connect Components Using MST Principles
  if(adjust_disconnection){
    while(length(filtered_components) > 1){
      cat(length(filtered_components), " Disconnected Components, apply MST.\n")
      
      # Find shortest edges between components
      edgesToAdd <- lapply(1:(length(filtered_components)-1), function(i) {
        lapply((i + 1):length(filtered_components), function(j) {
          findShortestEdge(filtered_components[[i]], filtered_components[[j]], feature_space)
        })
      })
      edgesToAdd = do.call(rbind.data.frame, lapply(unlist(edgesToAdd, recursive = FALSE),
                                                    function(x){data.frame(from = names(x$'from'),
                                                                           to = names(x$'to'),
                                                                           distance = x$distance)} )) %>% distinct()
      edges_group = edgesToAdd;
      edges_group[,1:2] = apply(edges_group[,1:2],2,function(y){
        apply(do.call(rbind,lapply(filtered_components,function(x) y %in% names(x))),2,function(x) which(x))
      })

      # Filter edges by MST
      g_meta <- igraph::graph_from_data_frame(edges_group, directed = FALSE)
      igraph::E(g_meta)$weight <- edges_group$distance

      mst_edges <- igraph::mst(g_meta, algorithm = "prim")
      mst_edges <- as.data.frame(igraph::get.edgelist(mst_edges))
      edgesToAdd_filtered = edgesToAdd[match(do.call(paste, mst_edges),do.call(paste, edges_group[,1:2])),]
      edgesToAdd = edgesToAdd_filtered
      
      # Add edges
      cat("add ", nrow(edgesToAdd), " edges\n")
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
    stop("Disconnected Components, Please increase knn or connect them.")
  }
  
  A = as.matrix(A)
  diag(A) = self_loop
  # Graph affinity matrix by symmetrize the adjacency mtx
  # W = (A + t(A))/2 # weighted
  W = pmin(A + t(A),1) # unweighted
  diag(W) = self_loop
  if(ncol(W) > 1e4){
    return(list(graph = as(W,"sparseMatrix"), adj_matrix = as(A,"sparseMatrix"), component = filtered_components))
  }else{
    return(list(graph = W, adj_matrix = A, component = filtered_components))
  }
}


#' Construct a Symmetric Gaussian Graph
#'
#' Constructs a Gaussian kernel graph from a given feature space.
#'
#' @param knn integer; the number of nearest neighbors to consider for each node. Default is 5.
#' @param feature_space matrix; a cell-by-coordinate matrix (e.g., 20 principal components)
#' @param alpha numeric; exponent used in the alpha-decaying Gaussian kernel
#' @param coef numeric; coefficient for adaptive bandwidth. Default is 1.
#' @param epsilon numeric; threshold below which edge weights are set to zero to sparsify the graph. Default is 1e-3.
#' @param self_loop numeric; value for the diagonal elements of the affinity matrix to add self-loops. Default is 1.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{graph}{A symmetric matrix representing the graph affinity matrix.}
#'   \item{adj_matrix}{A binary adjacency matrix indicating the presence of edges between nodes.}
#' }
#'
#' 
#' 
#' @keywords internal
#' 
#' @export
ConstructGaussianGraph <- function(knn = 5, feature_space, alpha = 10, coef = 1, 
                                   epsilon = 1e-3, self_loop = 1){
  cat("Constructing Gaussian kernel\n")
  node_num = nrow(feature_space)
  knn_list <- FNN::get.knn(feature_space, k = node_num - 1, algorithm = "kd_tree")
  Idx = knn_list$nn.index
  
  # Gaussian Kernel transfer
  knn_dist = knn_list$nn.dist
  bandwidth = coef * (knn_dist[,knn]) # adaptive bandwidth
  knn_dist = (knn_dist / bandwidth)^alpha
  knn_dist = exp(-knn_dist)
  knn_dist[knn_dist < epsilon] = 0 # make graph sparse
  
  A = Matrix::sparseMatrix(i = rep(1:node_num,each = ncol(Idx)),
                           j = c(t(Idx)),
                           dims = c(node_num, node_num),
                           x = c(t(knn_dist)))
  rownames(A) = colnames(A) = rownames(feature_space)
  A = as.matrix(A)
  diag(A) = self_loop
  
  # Graph affinity matrix (by default graph is connected)
  W = (A + t(A))/2
  A = W > 0
  
  return(list(graph = W,adj_matrix = A))
}


#' Construct a Symmetric SNN Graph
#'
#' Constructs a Shared Nearest Neighbor (SNN) graph from a given feature space using k-nearest neighbors (KNN) and Seurat's SNN computation.
#' 
#' @param knn integer; the number of nearest neighbors to consider for each node. Default is 5.
#' @param feature_space matrix; a cell-by-coordinate matrix (e.g., 20 principal components)
#' @param self_loop numeric; value for the diagonal elements of the affinity matrix to add self-loops. Default is 1.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{graph}{A symmetric matrix representing the graph affinity matrix.}
#'   \item{adj_matrix}{A binary adjacency matrix indicating the presence of edges between nodes.}
#' }
#'
#' 
#' 
#' @keywords internal
#' 
#' @export
ConstructSNNGraph <- function(knn = 5, feature_space, self_loop = 1, adjust_disconnection = TRUE){
  cat("Constructing SNN graph\n")
  knn_list <- FNN::get.knn(feature_space, k = knn, algorithm = "kd_tree")
  Idx = knn_list$nn.index
  Idx = cbind(1:nrow(Idx),Idx) # add self-loop
  rownames(Idx) = rownames(feature_space)
  snn_graph = Seurat:::ComputeSNN(Idx, prune = 0)
  snn_adj = as.matrix(snn_graph)
  rownames(snn_adj) <- colnames(snn_adj) <- rownames(feature_space)
  
  if(adjust_disconnection){
    snn_binary <- snn_adj > 0
    res = findDisconnectedComponents(snn_binary)
    
    # Remove singleton
    filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 10) comp else NULL)
    filtered_components <- Filter(Negate(is.null), filtered_components)
    filtered_node_names <- unlist(lapply(filtered_components, function(comp) names(comp)))
    snn_adj = snn_adj[rownames(snn_adj)%in%filtered_node_names,
          colnames(snn_adj)%in%filtered_node_names]
    if(dim(snn_adj)[1] < nrow(snn_binary)){
      cat("Remove ", nrow(snn_binary) - dim(snn_adj)[1], " singleton cells.\n")
    }
    
    while(length(filtered_components) > 1){
      cat(length(filtered_components), " Disconnected Components, apply MST.\n")
      
      # Find shortest edges between components
      edgesToAdd <- lapply(1:(length(filtered_components)-1), function(i) {
        lapply((i + 1):length(filtered_components), function(j) {
          findShortestEdge(filtered_components[[i]], filtered_components[[j]], feature_space)
        })
      })
      edgesToAdd = do.call(rbind.data.frame, lapply(unlist(edgesToAdd, recursive = FALSE),
                                                    function(x){data.frame(from = names(x$'from'),
                                                                           to = names(x$'to'),
                                                                           distance = x$distance)} )) %>% distinct()
      edges_group = edgesToAdd;
      edges_group[,1:2] = apply(edges_group[,1:2],2,function(y){
        apply(do.call(rbind,lapply(filtered_components,function(x) y %in% names(x))),2,function(x) which(x))
      })
      
      # Filter edges by MST
      g_meta <- igraph::graph_from_data_frame(edges_group, directed = FALSE)
      igraph::E(g_meta)$weight <- edges_group$distance
      
      mst_edges <- igraph::mst(g_meta, algorithm = "prim")
      mst_edges <- as.data.frame(igraph::get.edgelist(mst_edges))
      edgesToAdd_filtered = edgesToAdd[match(do.call(paste, mst_edges),do.call(paste, edges_group[,1:2])),]
      edgesToAdd = edgesToAdd_filtered
      
      # Add edges
      cat("add ", nrow(edgesToAdd), " edges\n")
      for (i in 1:nrow(edgesToAdd)) {
        idx_1 = which(rownames(Idx) == edgesToAdd[i,1])
        idx_2 = which(rownames(Idx) == edgesToAdd[i,2])
        node_ls = list(c(Idx[idx_1,],idx_2),
                       c(Idx[idx_2,],idx_1))
        weight = length(intersect(node_ls[[1]],node_ls[[2]])) / 
          length(union(node_ls[[1]],node_ls[[2]]))
        
        snn_adj[edgesToAdd[i,1],edgesToAdd[i,2]] = weight
        snn_adj[edgesToAdd[i,2],edgesToAdd[i,1]] = weight
      }
      snn_binary <- snn_adj > 0
      res = findDisconnectedComponents(snn_binary)
      filtered_components = res$components
    }
  }  
  
  # By definition, the adjacency mtx should be symmetric
  W = snn_adj
  diag(W) = self_loop
  A = W > 0
  if(ncol(W) > 1e4){
    return(list(graph = as(W,"sparseMatrix"), adj_matrix = as(A,"sparseMatrix")))
  }else{
    return(list(graph = W, adj_matrix = A))
  }
}

#' Construct Diffusion Operators
#'
#' Constructs a list of diffusion operators for a given symmetry affinity matrix of a graph.
#'
#' @param W matrix; symmetric affinity matrix of the graph.
#' @param max_time integer; the maximum diffusion time. The actual maximum diffusion time may be shorter if all nodes converge beforehand.
#'
#' @return A list of sparse matrices representing the diffusion operators at different times.
#' Each element in the list corresponds to P^t for t = 0, 2, 4, 8, ..., max_time.
#' The list is named by the corresponding diffusion time.
#'
#' @keywords internal
#' @export
ConstructDiffusionOperators <- function(W, max_time){
  cat("Create a list of diffusion operators...\n")
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  P = Doubly_stochastic(W)
  setTxtProgressBar(pb, 10)
  # P = RowwiseNormalize(W)
  # eig_res = RSpectra::eigs_sym(P, k = 1, which = "LM")
  max_step = max_time
  # P_ls = list(P)
  P_ls = NULL
  if(max_step < 1){
    stop("Incorrect diffusion time, no propogation")
  }else{
    t = 1
    # automatic decide max_step by checking whether the diag of P -> 1/n
    while(t <= floor(log(max_step,2)) & 
          max(abs((ncol(P) * diag(P)) * ncol(P) - ncol(P))) >= 1e-2 * ncol(P)){
      P = P %*% P; t = t + 1
      P_ls = c(P_ls,list(P))
      setTxtProgressBar(pb, 10 * (t + 1))
    }
    setTxtProgressBar(pb, 100)
  }
  # Add t = 0
  P_ls = c(list(diag(nrow(W))),P_ls)
  # make it sparse
  cat("\nConverting diffusion operators to sparse matrices...\n")
  P_ls = lapply(P_ls,function(x) as(x,"sparseMatrix"))
  names(P_ls) = c(0,2^seq(1,t-1))
  cat("\nMax diffusion time:",2^(t-1),"\n")
  return(P_ls)
}
