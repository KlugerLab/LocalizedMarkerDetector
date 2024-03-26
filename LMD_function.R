cran_packages <- c("igraph","ggplot2", "cowplot", "RColorBrewer", 
                   "data.table", "dplyr", "patchwork", 
                   "pheatmap", "ggplotify", "ggraph", 
                   "ClusterR", "Rcpp", "RcppArmadillo", 
                   "tictoc","svMisc","RSpectra","Matrix","progress","latex2exp","parallel")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)

# Mtx operation===========
# Enable C++11
Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# Create the C++ function
cppFunction('arma::mat fastMatMult(arma::mat A, arma::mat B) {
  arma::mat C = A * B;
  return C;
}', depends="RcppArmadillo")


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
  x = x[rowSums(x)!=0,,drop = FALSE]
  return( sweep(x, 1, rowSums(x), FUN = "/") )
}
is_sparse_matrix <- function(m) {
  any(grepl("Matrix", class(m)) & !grepl("dense", class(m)))
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
sinkhorn_knopp_largeData = function(A, niter = 100, tol = 1e-8, verb = FALSE) {
  # code refer to (https://rdrr.io/github/aaamini/nett/src/R/sinkhorn.R)
  delta = Inf
  for (irep in 1:niter) {
    scale_fac <- 1 / sqrt(rowSums(A))
    A <- Matrix::.sparseDiagonal(x = scale_fac) %*% A %*% Matrix::.sparseDiagonal(x = scale_fac)
    
    delta = pmax(max(abs(1 - rowSums(A))), max(abs(1 - colSums(A))))
    if (verb) nett::printf("err = %3.5e\n", delta)
    if (delta < tol){
      cat("large_graph, doubly-stochastic iter step: ",irep,"\n")
      break;
    }
  }
  return(A)
}

Doubly_stochastic <- function(W){
  if(ncol(W) < 1e4){
    scale_fac = sinkhorn_knopp(A = W, sym = TRUE)
    P = diag(scale_fac[[1]]) %*% W %*% diag(scale_fac[[2]])
    return(P)
  }else{
    return(sinkhorn_knopp_largeData(A = W))
  }
}

# Build Cell Graph & Transition Mtx =======
findDisconnectedComponents <- function(adj_matrix) {
  # Create a graph from the adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  # Find the connected components
  components <- igraph::components(g)
  # Get the number of disconnected components
  num_components <- components$no
  # Get the corresponding nodes for each component
  component_nodes <- split(V(g), components$membership)
  list(number_of_components = num_components, components = component_nodes)
}
findShortestEdge <- function(component1, component2, data) {
  # Calculate all pairwise distances between nodes in the two components
  distances <- outer(component1, component2, Vectorize(function(x, y) dist(data[c(x,y), ])))
  # Find the minimum distance and the corresponding nodes
  min_distance_idx <- which(distances == min(distances), arr.ind = TRUE)
  return(c(component1[min_distance_idx[1]], component2[min_distance_idx[2]]))
}

Symmetric_KNN_graph <- function(knn = 5, feature_space, adjust_by_MST = TRUE, self_loop = 1){
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
  filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 5) comp else NULL)
  filtered_components <- Filter(Negate(is.null), filtered_components)
  filtered_node_names <- unlist(lapply(filtered_components, function(comp) names(comp)))
  A = A[rownames(A)%in%filtered_node_names,
        colnames(A)%in%filtered_node_names]
  
  # Connect Components Using MST Principles
  if(adjust_by_MST){
    while(length(filtered_components) > 1){
      cat("Disconnected Components, adjust by MST.\n")
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
  if(ncol(W) > 1e4){
    return(list(graph = as(W,"sparseMatrix"), adj_matrix = as(A,"sparseMatrix"), component = filtered_components))
  }else{
    return(list(graph = W, adj_matrix = A, component = filtered_components))
  }
}

Symmetric_gaussian_graph <- function(knn = 5, feature_space, alpha = 1, coef = 2, epsilon = 1e-3, self_loop = 1){
  cat("Constructing Gaussian kernel\n")
  node_num = nrow(feature_space)
  knn_list <- FNN::get.knn(feature_space, k = node_num - 1, algorithm = "kd_tree")
  Idx = knn_list$nn.index

  # Gaussian Kernel transfer
  knn_dist = knn_list$nn.dist
  bandwidth = coef * (knn_dist[,knn])^2 # adaptive bandwidth
  knn_dist = (knn_dist^2 / bandwidth)^alpha
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

Obtain_Pls <- function(W, max_time){
  cat("Create a list of diffusion operator\n")
  pb <- txtProgressBar(min = 0, max = 100, style = 3)
  P = Doubly_stochastic(W)
  setTxtProgressBar(pb, 10)
  # P = Rowwise_normalize(W)
  # eig_res = RSpectra::eigs_sym(P, k = 1, which = "LM")
  max_step = max_time
  P_ls = list(P)
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
  names(P_ls) = 2^seq(0,t-1)
  # Add t = 0
  P_ls = c(list(diag(nrow(W))),P_ls) 
  names(P_ls)[1] = "0"
  cat("\nmax diffusion time:",2^(t-1))
  return(P_ls)
}
plot_knn_graph <- function(affinity_m, label = NULL, layout){
  layout <- as.matrix(layout)
  g <- graph_from_adjacency_matrix(affinity_m, 'undirected', weighted = TRUE)
  V(g)$frame.color <- NA
  
  if(!is.null(label)){
    V(g)$color <- label
    coldef <- setNames(
      colorRampPalette(brewer.pal(12, "Paired"))(length(unique(label))),
      unique(label))
    p = ggraph(g,layout = layout) + 
      geom_edge_link(aes(width = weight),color = "grey") + 
      geom_node_point(aes(color = color), size = 1) + 
      scale_edge_width_continuous(range = c(0.5, 1), breaks = c(0.5, 1)) + 
      scale_color_manual(values = coldef) +
      theme_graph()
  }else{
    p = ggraph(g,layout = layout) + 
      geom_edge_link(aes(width = weight),color = "grey") + 
      geom_node_point(size = 1) + 
      scale_edge_width_continuous(range = c(0.5, 1), breaks = c(0.5, 1)) + 
      theme_graph()
  }
  p = p + theme(
    plot.title = element_text(face="bold", size = 30),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15))
  return(p)
}

# Calculating Diffusion score =========
## Define a series of score_profile
get_score_profile <- function(state_0 = NULL, state_t = NULL, state_inf = NULL, P_diag_t = NULL,
                              score_ls = c("score0","delta_correction")){
  # score_ls = c("score0","score1","delta_correction","entropy")
  df = do.call(cbind,lapply(score_ls,function(s){
    if(s == "score0"){
      # KL(p^0||p^t)
      return(fastKLMatrix(state_0, state_t))
    }
    if(s == "max_score0"){
      return(fastKLVector(state_0,state_inf))
    }
    if(s == "score1"){
      # KL(p^t||p^inf)
      return(fastKLVector(state_t, state_inf))
    }
    if(s == "delta_correction"){
      # KL(p^0||p_delta^t)
      return(fastKLVector(state_0, P_diag_t))
    }
    if(s == "entropy"){
      # H(pt)
      return(-fastKLVector(state_t, rep(1,ncol(state_t))))
    }
  }))
  colnames(df) = score_ls
  return(df)
}

fast_calculate_score_profile <- function(W, max_time = 2^15, 
                                         init_state, P_ls = NULL,
                                         score_ls = c("score0")){
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  # final_state = rowSums(W)/sum(rowSums(W)) # row-stochastic
  final_state = rep(1/nrow(W),nrow(W)) # bi-stochastic

  # Calculate transition matrix
  if(is.null(P_ls)){
    P_ls = Obtain_Pls(W,max_time)
  }

  # Calculate multi-scale KL divergence
  score_df = do.call(cbind,lapply(1:length(P_ls),function(i){
    P = P_ls[[i]]
    state = fastMatMult(init_state, P)
    score_df = get_score_profile(state_0 = init_state,
                                 state_t = state,
                                 state_inf = final_state,
                                 # P_diag_t = diag(P),
                                 score_ls = score_ls)
    colnames(score_df) = paste0(colnames(score_df),"_",names(P_ls)[i])
    score_df
  }) )
  # Add scale_factor
  score_df = cbind(score_df,
                   get_score_profile(state_0 = init_state,
                                     state_t = init_state,
                                     state_inf = final_state,
                                     score_ls = c("entropy","max_score0")))
  rownames(score_df) = rownames(init_state)
  score_df = data.frame(score_df)
  return(score_df)
}
fast_calculate_score_profile_largeData <- function(W, max_time = 2^15, init_state,
                                                   score_ls = c("score0")){
  cat("Calculate LMD score profile for large data\n")
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  final_state = rep(1/nrow(W),nrow(W)) # bi-stochastic
  
  # Calculate transition matrix
  cat("Doubly Stochastic\n")
  P = Doubly_stochastic(W)
  # cat("Adjust self-loop weight")
  # diag(P) = 2/(min(rowSums(W)) - max(diag(W)))
  # P = Doubly_stochastic(P)
  
  # initialize score profile
  score_df = get_score_profile(state_0 = init_state,
                               state_t = init_state,
                               state_inf = final_state,
                               # P_t = diag(ncol(P)),
                               score_ls = score_ls)
  colnames(score_df) = paste0(colnames(score_df),"_","0")
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = floor(log(max_time,2)),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  # state_pre = init_state
  break_flag = FALSE
  for(t in 0:floor(log(max_time,2))){ # t here represents the exponent
    # cat("diffusion time:2^",t,"\n")
    # Automatic decide max_time by checking whether the diag of P -> 1/n
    if(!is_sparse_matrix(P) && (max(abs((ncol(P) * Matrix::diag(P)) * ncol(P) - ncol(P))) < 1e-2 * ncol(P))){
      cat("\nmax diffusion time:2^",t-1,"\n")
      break_flag = TRUE
      break
    }
    # Check the sparsity of P, if too dense, transfer it to dense matrix
    if(is_sparse_matrix(P) && (min(rowSums(P == 0)/ncol(P)) < 0.9)){
      P = as.matrix(P)
    }
    # cat("Graph dyadic\n")
    if(t > 0){
      if(is_sparse_matrix(P)){
        P = P %*% P
      }else{
        P = fastMatMult(P, P)
      }
    }
    state = fastMatMult(init_state, as.matrix(P))
    score_df_tmp = get_score_profile(state_0 = init_state,
                                 state_t = state,
                                 state_inf = final_state,
                                 # P_diag_t = diag(P),
                                 score_ls = score_ls)
    colnames(score_df_tmp) = paste0(colnames(score_df_tmp),"_",as.character(2^t))
    score_df = cbind(score_df,score_df_tmp)
    # score_df[,as.character(2^t)] = c(fastKLMatrix(init_state, state), fastKLVector(state, rep(1,ncol(state)))) # -entropy of pt
    pb$tick()
    # state_pre = state
  }
  if(!break_flag){cat("\nmax diffusion time:2^",t,"\n")}
  
  # Add scale_factor
  score_df = cbind(score_df,
                   get_score_profile(state_0 = init_state,
                                     state_t = init_state,
                                     state_inf = final_state,
                                     score_ls = c("entropy","max_score0")))
  
  rownames(score_df) = rownames(init_state)
  score_df = data.frame(score_df)
  return(score_df)
}
fast_calculate_score_profile_highres <- function(W, max_time = 100, init_state, 
                                                 score_ls = c("score0")){
  cat("Calculate LMD score profile for large data\n")
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  final_state = rep(1/nrow(W),nrow(W)) # bi-stochastic
  
  # Calculate transition matrix
  cat("Doubly Stochastic\n")
  P = Doubly_stochastic(W)
  P = as.matrix(P)
  score_df = get_score_profile(state_0 = init_state,
                               state_t = init_state,
                               state_inf = final_state,
                               # P_t = diag(ncol(P)),
                               score_ls = score_ls)
  colnames(score_df) = paste0(colnames(score_df),"_","0")
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = max_time,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  state = init_state
  break_flag = FALSE
  for(t in 1:max_time){
    # cat("diffusion time:2^",t,"\n")
    # Automatic decide max_time by checking whether the diag of P -> 1/n
    if(!is_sparse_matrix(P) && (max(abs((ncol(P) * Matrix::diag(P)) * ncol(P) - ncol(P))) < 1e-2 * ncol(P))){
      cat("\nmax diffusion time:",t-1,"\n")
      break_flag = TRUE
      break
    }
    state = fastMatMult(state, P)
    score_df_tmp = get_score_profile(state_0 = init_state,
                                     state_t = state,
                                     state_inf = final_state,
                                     # P_diag_t = diag(P),
                                     score_ls = score_ls)
    colnames(score_df_tmp) = paste0(colnames(score_df_tmp),"_",as.character(2^t))
    score_df = cbind(score_df,score_df_tmp)
    # score_df[,as.character(t)] = c(fastKLMatrix(init_state, state), fastKLVector(state, rep(1,ncol(state)))) # -entropy of pt
    pb$tick()
  }
  if(!break_flag){cat("\nmax diffusion time:",t,"\n")}
  
  # Add scale_factor
  score_df = cbind(score_df,
                   get_score_profile(state_0 = init_state,
                                     state_t = init_state,
                                     state_inf = final_state,
                                     score_ls = c("entropy","max_score0")))
  
  rownames(score_df) = rownames(init_state)
  score_df = data.frame(score_df)
  
  return(score_df)
}
Calculate_outgoing_weight <- function(W, max_time = 2^15, init_state){
  cat("Calculate LMD score profile for large data\n")
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  final_state = rep(1/nrow(W),nrow(W)) # bi-stochastic
  
  # Calculate transition matrix
  cat("Doubly Stochastic\n")
  P = Doubly_stochastic(W)
  cat("self-weight:",diag(P)[1],"\n")
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = floor(log(max_time,2)),
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  indicate_mtx = init_state > 0
  all_weight_to_detect = t(W %*% t(indicate_mtx))
  out_weight = rowSums(all_weight_to_detect * (!indicate_mtx))
  inner_weight = (rowSums(all_weight_to_detect * indicate_mtx) + rowSums(indicate_mtx))/2
  out_in_ratio = out_weight / inner_weight
  out_in_ratio_df = data.frame("0" = c(inner_weight,rowSums(indicate_mtx)))
  break_flag = FALSE
  for(t in 0:floor(log(max_time,2))){
    # cat("diffusion time:2^",t,"\n")
    # Automatic decide max_time by checking whether the diag of P -> 1/n
    if(!is_sparse_matrix(P) && (max(abs((ncol(P) * Matrix::diag(P)) * ncol(P) - ncol(P))) < 1e-2 * ncol(P))){
      cat("\nmax diffusion time:2^",t-1,"\n")
      break_flag = TRUE
      break
    }
    # Check the sparsity of P, if too dense, transfer it to dense matrix
    if(is_sparse_matrix(P) && (min(rowSums(P == 0)/ncol(P)) < 0.9)){
      P = as.matrix(P)
    }
    # cat("Graph dyadic\n")
    if(t > 0){
      if(is_sparse_matrix(P)){
        P = P %*% P
      }else{
        P = fastMatMult(P, P)
      }
    }
    state = fastMatMult(init_state, as.matrix(P))
    pb$tick()
    
    # Add a checkpoint of outgoing weight
    indicate_mtx = state > 0
    all_weight_to_detect = t(W %*% t(indicate_mtx))
    out_weight = rowSums(all_weight_to_detect * (!indicate_mtx))
    inner_weight = (rowSums(all_weight_to_detect * indicate_mtx) + rowSums(indicate_mtx))/2
    out_in_ratio = out_weight / inner_weight
    out_in_ratio_df[,as.character(2^t)] = c(inner_weight,rowSums(indicate_mtx))
  }
  colnames(out_in_ratio_df)[1] = "0"
  if(!break_flag){cat("\nmax diffusion time:2^",t-1,"\n")}
  out_in_ratio_df = cbind(out_in_ratio_df[1:nrow(init_state),],
                          out_in_ratio_df[(nrow(init_state)+1):nrow(out_in_ratio_df),])
  colnames(out_in_ratio_df) = paste0(rep(c("ratio_","count_"), each = ncol(out_in_ratio_df)/2),
                              colnames(out_in_ratio_df))
  ## Entropy
  out_in_ratio_df[,"entropy"] = -fastKLVector(init_state, rep(1,ncol(init_state)))
  rownames(out_in_ratio_df) = rownames(init_state)
  return(out_in_ratio_df)
}

obtain_lmds = function(score_profile, correction = FALSE){
  # get LMDS
  sub_score0 = grep("^score0",colnames(score_profile))
  sub_delta_correction = grep("^delta_correction",colnames(score_profile))
  if(!correction){
    # raw LMDS
    # Scale the score profile by maximum value
    df = score_profile/score_profile[,"max_score0"]
    cumulative_score = rowSums(df[,sub_score0])
  }else{
    # adjusted LMDS
    # Scale the score profile by entropy
    df = score_profile/score_profile[,"entropy"]
    cumulative_score = rowSums(df[,sub_score0] - df[,sub_delta_correction])
  }
  return(cumulative_score)
}

# combine fast_calculate_score_profile & obtain_lmds in one step
fast_get_lmds <- function(W, max_time = 2^15, init_state, P_ls = NULL, correction = FALSE, largeData = TRUE, highres = FALSE){
  if(highres){
    score_profile = fast_calculate_score_profile_highres(W = W, max_time = max_time, 
                                                           init_state = init_state)
    cumulative_score = obtain_lmds(score_profile = score_profile, correction = correction)
    return(list(score_profile = score_profile,cumulative_score = cumulative_score))
  }
  if((ncol(W) > 1e4) | largeData){
    score_profile = fast_calculate_score_profile_largeData(W = W, max_time = max_time, 
                                                 init_state = init_state)
  }else{
    score_profile = fast_calculate_score_profile(W = W, max_time = max_time, 
                                                 init_state = init_state, P_ls = P_ls)
  }
  cumulative_score = obtain_lmds(score_profile = score_profile, correction = correction)
  return(list(score_profile = score_profile,cumulative_score = cumulative_score))
}

knee_point = function(Vecs){
  Vecs = sort(Vecs)
  # calculating the distance of each point on the curve from a line drawn from the first to the last point of the curve. The point with the maximum distance is typically considered the knee point.
  curve_data = as.matrix(data.frame(x = 1:length(Vecs), y = Vecs))
  line_start <- curve_data[1, ]
  line_end <- curve_data[nrow(curve_data), ]
  
  dx <- line_end["x"] - line_start["x"]
  dy <- line_end["y"] - line_start["y"]
  norm_factor <- sqrt(dx^2 + dy^2)
  A <- dy
  B <- -dx
  C <- dx * line_start["y"] - dy * line_start["x"]
  distances <- abs(A * curve_data[, "x"] + B * curve_data[, "y"] + C) / norm_factor
  knee_point = which.max(distances)
  
  plot(Vecs,pch = 20)
  points(knee_point,Vecs[knee_point],col='red',pch = 20)
  return(knee_point)
}
LMD <- function(expression, feature_space, knn = 5, 
                kernel = FALSE, max_time = 2^15, adjust_bridge = TRUE, self_loop = 1,
                score_correction = FALSE, largeData = TRUE, highres = FALSE, min_cell = 5){
  if(any(colnames(expression) != rownames(feature_space))){stop("Cells in expression mtx and feature space don't match.")}
  if(kernel){
    W = Symmetric_gaussian_graph(knn = knn, feature_space = feature_space, alpha = 1, coef = 2, epsilon = 1e-3, self_loop = self_loop)$'graph'
  }else{
    W = Symmetric_KNN_graph(knn = knn, feature_space = feature_space, adjust_by_MST = adjust_bridge, self_loop = self_loop)$'graph'
  }
  rho = Rowwise_normalize(expression)
  rho = rho[which(apply(rho,1,function(x) sum(x>0) >= min_cell))
            ,,drop = FALSE] # sanity check & remove genes which express at less than 5 cells
  cat(sprintf("Remove %d genes which express at less than %d cells\n",
              nrow(expression) - nrow(rho),min_cell))
  res = fast_get_lmds(W = W, max_time = max_time,
                                   init_state = rho, 
                      correction = score_correction, 
                      largeData = largeData,
                      highres = highres)
  return(res)
}
show_result_lmd <- function(res.lmd, n = length(res.lmd$cumulative_score)){
  score = res.lmd$cumulative_score
  score = sort(score)
  df = data.frame(score = score)
  df$'rank' = 1:nrow(df)
  # print(head(df,n = n))
  gene_rank = setNames(df$'rank',rownames(df))
  return(list(gene_table = df, gene_rank = gene_rank, cut_off_gene = gene_rank[1:knee_point(score)]))
}

# Visualization =========
FeaturePlot_custom <- function(coord, value, reduction = NULL, dims = 1:2, value_name = NULL,title_name = NULL,order_point = TRUE){
  if(class(coord) == "Seurat"){
    if(is.character(value)){
      if(value %in% colnames(coord@meta.data)){
        value = coord@meta.data[,value]
      }
      else{
        stop("Feature is not found in the given seurat object")
      }
    }
    
    if(!is.null(reduction)){
      coord <- Embeddings(coord, reduction = reduction)[,dims]
    } else if ("umap" %in% names(coord@reductions)){
      coord <- Embeddings(coord, reduction = "umap")[,dims]
    } else if ("tsne" %in% names(coord@reductions)){
      coord <- Embeddings(coord, reduction = "tsne")[,dims]
    } else{
      stop("Neither UMAP nor t-SNE embeddings are found in the Seurat object.")
    }
  }else if(is.null(coord)){
    stop("coord is missing")
  }
  if(length(value)!=nrow(coord)){stop("Unmatched Dimension!")}
  df = data.frame(cbind(coord[,1:2],value = value))
  if(order_point){df = df[order(df[,3]),]}
  p <- ggplot(df, aes(x=!!as.name(colnames(coord)[1]), y=!!as.name(colnames(coord)[2]), color=value)) + 
    geom_point(size=0.5)
  if(!is.factor(df[,3])){
    p <- p + scale_color_gradient(low = "lightgrey", high = "blue")
  }
  p = p + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
    theme(
      plot.title = element_text(face="bold", size = 30),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 15),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank()) + 
    labs(color = value_name, title = title_name)
  return(p)
  
}
FeaturePlot_diffusion <- function(coord, init_state, P_ls = NULL, W = NULL, check_time = NULL, gene_name = NULL, gene_color = "blue"){
  init_state = as.matrix(init_state)
  tmp <- if (!is.null(P_ls) && length(P_ls) > 0) {
    P_ls[[1]]
  } else if (!is.null(W)) {
    W
  } else {
    NULL
  }
  if(is.null(tmp) & is.null(W)){stop("Cell Graph 'W' not found")}
  graph_node_number = ncol(tmp); rm(tmp)
  if((ncol(init_state)!=graph_node_number) & (nrow(init_state) == graph_node_number)){
    init_state = t(init_state)
  }
  if(ncol(init_state)!=graph_node_number | nrow(init_state)!=1){stop("Unmatched dimension, can only take expression vector for one gene!")}
  check_step = sort(unique(floor(log(check_time,base = 2))))
  check_step = check_step[check_step >= 0]
  
  if(0 %in% check_time){
    check_time = c(0,2^check_step)
  }else{
    check_time = 2^check_step
  }
  if(is.null(P_ls)){P_ls = Obtain_Pls(W,max(check_time))}
  check_time = check_time[check_time %in% names(P_ls)]
  sub = match(check_time,names(P_ls))
  multi_state = sapply(P_ls[sub],function(P){
    state = fastMatMult(init_state, P)
  })
  
  colnames(multi_state) = check_time
  # degree_node = rowSums(W) # Normalize multi_state with degree of node for visualization
  pl = lapply(seq(ncol(multi_state)),function(i){
    legend_break_label = seq(0, max(multi_state[,i]), length.out = 2)
    if(max(legend_break_label)>1e-3){
      legend_break_label = signif(legend_break_label, digits = 2)
    }else{
      legend_break_label = formatC(legend_break_label, format = "e", digits = 0)
    }
    p = suppressWarnings({
      FeaturePlot_custom(coord = coord, value = multi_state[,i]/max(multi_state[,i]),
                                 title_name = ifelse(as.numeric(colnames(multi_state)[i])==0,
                                                     latex2exp::TeX("$T = 0$"),
                                                     latex2exp::TeX(paste0("$T = 2^{",log(as.numeric(colnames(multi_state)[i]),base = 2),"}$"))),
                                 order_point = TRUE) + 
        scale_color_gradient(name = "Density",
                             low = "lightgrey", high = gene_color,
                             limits = c(0,1),
                             breaks = seq(0, 1, length.out = 2),
                             labels = legend_break_label)
    })
    p
  })
  
  p = wrap_plots(pl, nrow = 1) + plot_annotation(title = gene_name,
                                                 theme = theme(plot.title = element_text(face="bold",size = 25)) )
  
  return(p)
}
FeaturePlot_meta <- function(dat, coord = NULL, feature_partition, reduction = NULL, assays = "RNA"){
  if(class(dat) == "Seurat"){
    if (!is.null(reduction)){
      coord <- Embeddings(dat, reduction = reduction)
    } else if ("umap" %in% names(dat@reductions)){
      coord <- Embeddings(dat, reduction = "umap")
    } else if ("tsne" %in% names(dat@reductions)){
      coord <- Embeddings(dat, reduction = "tsne")
    } else{
      stop("Neither UMAP nor t-SNE embeddings are found in the Seurat object.")
    }
    dat = as.matrix(dat[[assays]]@data)
  }else if(is.null(dat)){
    stop("Data missing")
  }else if(is.null(coord)){
    stop("coordinate missing")
  }
  feature_partition = as.factor(feature_partition)
  pl <- lapply(levels(feature_partition), function(level){
    genes = names(feature_partition)[feature_partition == level]
    plot.title1 <- sprintf("Module %s (%d genes)",level,length(genes))
    df = data.frame(cbind(coord[colnames(dat),,drop = FALSE],
               value = apply(dat[genes,,drop = FALSE],2,mean)))
    p1 <- ggplot(df[order(df$value),], aes(x=!!as.name(colnames(df)[1]), y=!!as.name(colnames(df)[2]), color=value)) + 
      geom_point(size = 0.2) + 
      scale_color_gradient(low = "lightgrey", high = "blue") + 
      ggtitle(plot.title1) + labs(color = "Expression") +
      theme(legend.title = element_text(size = 8),
            panel.grid = element_blank(),
            panel.background = element_blank()) 
    p1
  })
  names(pl) = levels(feature_partition)
  return(pl)
}
Visualize_score_pattern <- function(score_profile, genes = NULL, 
                                    label_class = NULL, facet_class = NULL, 
                                    add_point = NULL, dyadic = TRUE, text = FALSE){
  score_df = score_profile
  if(!all(genes %in% rownames(score_df))){stop("Genes not found!")}
  profiles = names(which(table(sub("_\\d+", "", colnames(score_df)))>1))
  score_df = score_df[genes,]
  score_df$'gene' = rownames(score_df)
  score_df$'label' = score_df$'gene'
  score_df$'facet' = NA
  
  if(!is.null(label_class)){
    if(length(label_class)!=length(genes)){
      warning("Labels doesn't match genes, ignore labels")
    }else{
      score_df$'label' = label_class
    }
  }
  
  if(!is.null(facet_class)){
    if(length(facet_class)!=length(genes) & facet_class!="profiles"){
      warning("Facet labels doesn't match genes, ignore facet labels")
    }else{
      score_df$'facet' = facet_class
    }
  }

  df <- reshape2::melt(score_df,
                       id = colnames(score_df)[!grepl(paste0(profiles,"_",collapse = "|"),colnames(score_df))],
                       variable.name = c("step"), value.name = "score")
  df$profiles = sub("_([^_]*)$", "\\1", sub("\\d+$", "", df$step))
  df$step = sub(".*_([0-9]+).*", "\\1", df$step)
  x_breaks = unique(as.numeric(df$step))
  if(max(diff(x_breaks)) > 1 & dyadic){
    names(x_breaks) = log(x_breaks,base = 2)
    names(x_breaks) = ifelse(names(x_breaks) == "-Inf","0",paste0("2^",names(x_breaks)))
    x_breaks = setNames(names(x_breaks),x_breaks)
    df$step = as.factor(as.numeric(df$step))
    levels(df$step) = x_breaks[levels(df$step)]
  }else{
    df$step = as.numeric(df$step)
  }
  if(!is.null(label_class)){
    p = ggplot(data = df, mapping = aes(x = step, y = score, color = label, linetype = profiles)) + 
      geom_line(aes(group = interaction(gene,profiles))) + labs(x = "Time", y = "Normalized Diffusion KL Score")
    if(text == TRUE){
      p = p +
        geom_text(data = df %>% group_by(interaction(gene,profiles)) %>% slice_head(n = 4) %>% slice_tail(n = 1),
                  aes(label = gene), hjust = 0) + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }

  }else{
    p = ggplot(data = df, mapping = aes(x = step, y = score, color = gene, linetype = profiles)) + 
      geom_line(aes(group = interaction(gene,profiles))) + labs(x = "Time", y = "Normalized Diffusion KL Score")
  }
  if(!is.null(add_point)){
    add_point = x_breaks[as.character(add_point)]
    add_point = add_point[!is.na(add_point)]
    p = p + geom_point(data = df %>% filter(step %in% add_point), aes(x = step, y = score, color = label), size = 2) + 
      theme(axis.text.x = element_text(
        angle = 45,hjust = 1,
        color = ifelse(levels(df$step) %in% add_point, "red", "black"),
        face = ifelse(levels(df$step) %in% add_point, "bold", "plain")
      ))
  }
  if(!is.null(facet_class)){
    if(facet_class == "profiles"){
      p = p + facet_wrap(~profiles,scales = "free")
    }else{
      p = p + facet_wrap(~facet,scales = "free")
    }
  }
  p = p + theme(
    plot.title = element_text(face="bold", size = 30),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    # panel.grid = element_blank(),
    # panel.background = element_blank(),
    axis.line = element_line(colour = "black")
    ) 
  # + ylim(0,1)
  return(p)
}
Visualize_jaccard_mtx <- function(jaccard_mtx){
  if(is.null(rownames(jaccard_mtx))|is.null(colnames(jaccard_mtx))){
    stop("Please define celltype names & module names")
  }
  df = reshape2::melt(jaccard_mtx)
  colnames(df)[1:2] = c("CellType","Module")
  df$CellType = factor(df$CellType,levels = rownames(jaccard_mtx))
  df$Module = factor(df$Module,levels = rev(colnames(jaccard_mtx)))
  
  p <- ggplot(df, aes(x=CellType, y=Module, fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradientn(colors=colorRampPalette(rev(c("firebrick3", "white")))(99),limits = c(0,1), name="Jaccard Index") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1,size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 15),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + geom_text(data = df %>% filter(value > 0.4), 
                  aes(label = sprintf("%.2f", value)), 
                  color = "white", size = 4)
  return(p)
}

# Obtain Gene Modules ========
## Calculate gene pairwise distance
Calculate_distance <- function(dat, method){
  rho_dat = Rowwise_normalize(dat)
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

# Obtain Modules Activity Score ========
## Calculate Cell Module-activity score
Obtain_cell_partition <- function(expr_dat, gene_partition, 
                                  cell_kNN_graph = NULL, major_vote = 5){
  # Fit GMM
  if(sum(!names(gene_partition) %in% rownames(expr_dat))){
    stop("gene name doesn't match data")}
  
  meta_g = sapply(levels(gene_partition),function(topic){
    colMeans(expr_dat[names(gene_partition)[gene_partition == topic],,drop = FALSE])
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
  labels = predict(gmm,as.matrix(vector))
  
  
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
FindPC = function(srat){
  stdv <- srat[["pca"]]@stdev
  sum.stdv <- sum(srat[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  return(min.pc)
}
AddModuleActivityScore <- function(srat, gene_partition, assay = "RNA", do_local_smooth = TRUE, knn = 10, major_vote = 5, nloop = 100){
  # adjust the format of gene_partition
  gene_partition = setNames(as.character(gene_partition),names(gene_partition))
  gene_partition = gene_partition[names(gene_partition) %in% rownames(srat)]
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
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
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
  colnames(cell_block) = paste0("Module",colnames(cell_block))
  srat <- AddMetaData(srat, cell_block, col.name = colnames(cell_block))
  return(srat)
}


# Test Function =============
Generate_Pseudo_gene <- function(data){
  randseed = 233
  set.seed(randseed)
  seed <- sample(c(1:1e5),size=1e5)
  
  dat.detection = data > apply(data,2,median)
  each.gene.detection = apply(dat.detection,1,sum)
  
  # sample genes based on the hist of each.gene.detection
  xhist=hist(each.gene.detection,breaks = 1000,plot = FALSE)
  each.gene.bin = ceiling(each.gene.detection / unique(diff(xhist$breaks)))
  
  # choose bins where each sample belongs to
  samplesize=10000
  set.seed(seed[1])
  bins=with(xhist,sample(length(mids),samplesize,p=density,replace=TRUE))
  
  # sample genes based on the bins
  gene_sub = unlist(
    lapply(unique(bins),function(x, i){
      set.seed(seed[i])
      sample_num = table(bins)[which(names(table(bins)) == x)]
      sample(names(each.gene.bin)[each.gene.bin == x])[1:sample_num]
    }, i = 1:length(unique(bins)))
  )
  gene_sub = gene_sub[!is.na(gene_sub)]
  
  mat = data[gene_sub,]
  pseudo_dat = do.call(rbind,purrr::map2(seq_len(nrow(mat)), seed[1:nrow(mat)], ~ {
    set.seed(.y)
    sample(mat[.x, ])
  }) )
  
  colnames(pseudo_dat) = colnames(data)
  rownames(pseudo_dat) = paste0("pseudo-",gene_sub)
  
  return(pseudo_dat)
}

quick_marker_benchmark <- function(gene_rank_vec,
                                   folder_path = "/data/ruiqi/local_marker/LocalMarkerDetector/benchmark_result",
                                   tissue_name = "marrow_facs"){
  max_logfc = read.table(file.path(folder_path, paste0(tissue_name,"_ground_truth_c1.txt"))) %>% rownames()
  celldb_marker = read.table(file.path(folder_path, paste0(tissue_name,"_ground_truth_c2.txt")))[,1]
  gt_list = c(lapply(seq(50,1000,50),function(x) max_logfc[1:x]),list(celldb_marker))
  names(gt_list) = c(paste0("Top",seq(50,1000,50)),"CellMarkerDB")
  df_benchmark = data.frame(gene_rank_vec,row.names = names(gene_rank_vec))
  
  auc_vec = do.call(c,lapply(1:length(gt_list),function(i){
    true_marker = gt_list[[i]]
    df_benchmark$"gt" = 0
    df_benchmark[true_marker,"gt"] = 1
    
    library(pROC)
    roc = roc(df_benchmark$gt, df_benchmark[,1], direction = ">")
    as.numeric(auc(roc))
  }))
  names(auc_vec) = names(gt_list)
  return(auc_vec)
}


