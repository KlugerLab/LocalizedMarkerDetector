library("Rcpp")
library("RcppArmadillo")
library("igraph")

# FUNCTION ====================
# Mtx operation
# Enable C++11
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# Create the C++ function
## Fast Matrix Multiplication
cppFunction('arma::mat fastMatMult(arma::mat A, arma::mat B) {
  arma::mat C = A * B;
  return C;
}', depends="RcppArmadillo")

## Fast derivation of KL between pair-wise columns of two matrix (eg. KL between 1st col from mtxA and 1st col from mtxB)
Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;

NumericVector calcKL(NumericVector p, NumericVector q) {
  double sum = 0;
  for (int j = 0; j < p.size(); ++j) {
    if(p[j] == 0) continue;
    double val = p[j] * log(p[j] / q[j]);
    if (NumericVector::is_na(val) || std::isinf(val)) val = 0.0;
    sum += val;
  }
  return NumericVector::create(sum);
}

// [[Rcpp::export]]
NumericVector fastKLMatrix(NumericMatrix x, NumericMatrix y) {
  int n = x.nrow();
  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    out[i] = calcKL(x.row(i), y.row(i))[0];
  }
  return out;
}

// [[Rcpp::export]]
NumericVector fastKLVector(NumericMatrix x, NumericVector y) {
  int n = x.nrow();
  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    out[i] = calcKL(x.row(i), y)[0];
  }
  return out;
}
')
Rowwise_normalize <- function(x){
  return( sweep(x, 1, rowSums(x), FUN = "/") )
}

# Doubly Stochastic
l2_norm = function(x) sqrt(sum(x^2))
sinkhorn_knopp = function(A, sums = rep(1, nrow(A)),
                          niter = 100, tol = 1e-8, sym = FALSE, verb = FALSE) {
  # code refer to (https://rdrr.io/github/aaamini/nett/src/R/sinkhorn.R)
  delta = Inf
  r = c = rep(1, nrow(A))
  converged = FALSE
  t = 1
  while( t <= niter && !converged) {
    r = sums / (A %*% c)
    cnew = sums / (t(A) %*% r)
    
    # Symmetric Sinkhorn--Knopp algorithm could oscillate between two sequences,
    # need to bring the two sequences together (See for example "Algorithms For
    # The Equilibration Of Matrices And Their Application To Limited-memory
    # Quasi-newton Methods")
    if (sym) cnew = (cnew + r)/2
    
    delta = l2_norm(cnew-c)
    if (delta < tol) converged = TRUE
    if (verb) nett::printf("err = %3.5e\n", delta)
    c = cnew
    t = t+1
  }
  list(r = as.numeric(r), c = as.numeric(c))
}

Doubly_stochastic <- function(W){
  scale_fac = sinkhorn_knopp(A = W, sym = TRUE)
  P = diag(scale_fac[[1]]) %*% W %*% diag(scale_fac[[2]])
  
  return(P)
}

# Build Cell Graph
## Check disconnected components in graph
findDisconnectedComponents <- function(adj_matrix) {
  # Create a graph from the adjacency matrix
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  # Find the connected components
  components <- components(g)
  # Get the number of disconnected components
  num_components <- components$no
  # Get the corresponding nodes for each component
  component_nodes <- split(V(g), components$membership)
  list(number_of_components = num_components, components = component_nodes)
}
## Connect two disconnected components by shortest edge
findShortestEdge <- function(component1, component2, data) {
  # Calculate all pairwise distances between nodes in the two components
  distances <- outer(component1, component2, Vectorize(function(x, y) dist(data[c(x,y), ])))
  # Find the minimum distance and the corresponding nodes
  min_distance_idx <- which(distances == min(distances), arr.ind = TRUE)
  return(c(component1[min_distance_idx[1]], component2[min_distance_idx[2]]))
}
## Construct kNN graph
Symmetric_KNN_graph <- function(knn = 5, feature_space, adjust_by_MST = TRUE){
  knn_list <- FNN::get.knn(feature_space, k = knn, algorithm = "kd_tree")
  Idx = knn_list$nn.index
  A = Matrix::sparseMatrix(i = rep(1:nrow(Idx),each = knn),
                           j = c(t(Idx)),
                           dims = c(nrow(Idx), nrow(Idx)),
                           x = 1)
  rownames(A) = colnames(A) = rownames(feature_space)
  
  # remove orphans/singletons
  res = findDisconnectedComponents(A)
  filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 5) comp else NULL)
  filtered_components <- Filter(Negate(is.null), filtered_components)
  filtered_node_names <- unlist(lapply(filtered_components, function(comp) names(comp)))
  A = A[rownames(A)%in%filtered_node_names,
        colnames(A)%in%filtered_node_names]
  
  # Connect Components Using MST Principles
  if(adjust_by_MST){
    while(length(filtered_components) > 1){
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
  
  if (length(filtered_components) >= 2) {
    stop("Disconnected Components, Please increase knn.")
  }
  
  A = as.matrix(A)
  
  # Graph affinity matrix by symmetrize the adjacency mtx
  # W = (A + t(A))/2 # weighted
  W = pmin(A + t(A),1) # unweighted
  
  return(list(graph = W,adj_matrix = A, component = filtered_components))
}
Symmetric_gaussian_graph <- function(knn = 5, feature_space, alpha = 1, coef = 2, epsilon = 1e-3){
  node_num = nrow(feature_space)
  knn_list <- FNN::get.knn(feature_space, k = node_num - 1, algorithm = "kd_tree")
  Idx = knn_list$nn.index
  
  # Gaussian Kernel transfer
  knn_dist = knn_list$nn.dist
  bandwidth = coef * (knn_dist[,knn])^2
  knn_dist = (knn_dist^2 / bandwidth)^alpha
  knn_dist = exp(-knn_dist)
  knn_dist[knn_dist < epsilon] = 0 # make graph sparse
  
  A = Matrix::sparseMatrix(i = rep(1:node_num,each = ncol(Idx)),
                           j = c(t(Idx)),
                           dims = c(node_num, node_num),
                           x = c(t(knn_dist)))
  rownames(A) = colnames(A) = rownames(feature_space)
  A = as.matrix(A)
  
  # Graph affinity matrix (by default graph is connected)
  W = (A + t(A))/2
  A = W > 0
  
  return(list(graph = W,adj_matrix = A))
}

# Save a series of P^t (t = 2^1, 2^2, 2^3, ...) as a large list
Obtain_Pls <- function(W, max_time){
  P = Doubly_stochastic(W)
  # P = Rowwise_normalize(W)
  # eig_res = RSpectra::eigs_sym(P, k = 1, which = "LM")
  max_step = max_time
  P_ls = NULL
  if(max_step < 1){
    print("no propogation")
  }else{
    t = 1
    while(t <= floor(log(max_step,2))){
      P = P %*% P; t = t + 1
      P_ls = c(P_ls,list(P))
    }
  }
  names(P_ls) = 2^seq(1,floor(log(max_step,2)))
  return(P_ls)
}

# Calculate KL-score & pulse-score & LMDS score
fast_calculate_multi_score <- function(W, max_time = 2^15, init_state, P_ls = NULL, correction = TRUE){
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  # degree = rowSums(W)/sum(rowSums(W)) # row-stochastic
  degree = rep(1/nrow(W),nrow(W)) # bi-stochastic
  
  # Calculate transition matrix
  if(is.null(P_ls)){
    P_ls = Obtain_Pls(W,max_time)
  }
  P_ls = c(list(diag(nrow(W))),P_ls) # Add t = 0
  names(P_ls)[1] = "0"
  
  # Calculate multi-scale KL divergence
  score_df = do.call(cbind,lapply(P_ls,function(P){
    state = fastMatMult(init_state, P)
    c(fastKLMatrix(init_state, state), fastKLVector(init_state, diag(P)))
  }) )
  score_df = cbind(score_df[1:nrow(init_state),],
                   score_df[(nrow(init_state)+1):nrow(score_df),])
  colnames(score_df) = paste0(rep(c("score0_","correction_"), each = ncol(score_df)/2),
                              colnames(score_df))
  # Maximum value
  score_df = cbind(score_df,fastKLVector(init_state,degree))
  colnames(score_df)[ncol(score_df)] = "maximum_val"
  rownames(score_df) = rownames(init_state)
  
  # Normalized the score profile
  score_df = data.frame(score_df/score_df[,"maximum_val"])
  
  # get LMDS
  sub_score0 = grep("score0",colnames(score_df))
  sub_correction = grep("correction",colnames(score_df))
  if(correction){
    cumulative_score = rowSums(score_df[,sub_score0] - score_df[,sub_correction])
  }else{
    cumulative_score = rowSums(score_df[,sub_score0])
  }
  
  return(list(score_profile = score_df,cumulative_score = cumulative_score))
}

# Get LMDS in one step
LMD <- function(expression, feature_space, knn = 5, 
                kernel = FALSE, max_time = 2^15, adjust_bridge = TRUE,
                score_correction = TRUE){
  if(any(colnames(expression) != rownames(feature_space))){stop("Cells in expression mtx and feature space don't match.")}
  if(kernel){
    W = Symmetric_gaussian_graph(knn = knn, feature_space = feature_space, alpha = 1, coef = 2, epsilon = 1e-3)$'graph'
  }else{
    W = Symmetric_KNN_graph(knn = knn, feature_space = feature_space, adjust_by_MST = adjust_bridge)$'graph'
  }
  rho = Rowwise_normalize(expression)
  res = fast_calculate_multi_score(W = W, max_time = max_time,
                                   init_state = rho, correction = score_correction)
  return(res)
}
show_result_lmd <- function(res.lmd, n = length(res.lmd$cumulative_score)){
  score = res.lmd$cumulative_score
  score = sort(score)
  df = data.frame(score = score)
  df$'rank' = 1:nrow(df)
  print(head(df,n = n))
  gene_rank = setNames(df$'rank',rownames(df))
  return(list(gene_table = df, gene_rank = gene_rank, cut_off_gene = gene_rank[1:knee_point(score)]))
}


# Step for calculating LMDS ====================
#' feature_space: 20PC coordinates
#' dat: bone marrow log-transformed count data

#' Build KNN Graph
knn_result = Symmetric_KNN_graph(knn = 5, feature_space = feature_space)
A = knn_result$adj_matrix # Adjacency Matrix
W = knn_result$graph # Symmetrized Graph ((A + AT) / 2)

#' Generate Generate $P^{t}, t = 2,4,8,...,max_time$
P_ls = Obtain_Pls(W = W, max_time = 2^15)

#' Obtain score-profile & pulse-score profile & LMDS
score_result = fast_calculate_multi_score(W = W, init_state = Rowwise_normalize(dat), P_ls = P_ls, 
                           correction = TRUE) # If correction TRUE, return adjusted LMDS; else return LMDS


