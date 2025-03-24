# =====================================
# Calculating LMD Score ===============
# =====================================
# Extrinsic Function: fast_get_lmds, LMD, show_result_lmd

#' Compute Score Profile
#'
#' Computes various scores for state transitions based on different metrics such as KL divergence and entropy.
#'
#' @param state_0 matrix; the initial state matrix.
#' @param state_t matrix; the state matrix at time t.
#' @param state_inf matrix; the state matrix at equilibrium.
#' @param P_diag_t matrix; the diagonal elements of the diffusion operator at time t.
#' @param score_ls character vector; a list of scores to compute. Possible values are "score0", "max_score0", "score1", "delta_correction", and "entropy". Default is \code{c("score0", "delta_correction")}.
#' \itemize{
#'   \item "score0": KL divergence between the initial state (\code{state_0}) and the state at time t (\code{state_t}).
#'   \item "max_score0": KL divergence between the initial state (\code{state_0}) and the equilibrium state (\code{state_inf}).
#'   \item "score1": KL divergence between the state at time t (\code{state_t}) and the equilibrium state (\code{state_inf}).
#'   \item "delta_correction": KL divergence between the initial state (\code{state_0}) and the diagonal elements of the diffusion operator at time t (\code{diag(P^t)}).
#'   \item "entropy": Entropy of the state at time t (\code{state_t}).
#' }
#' @return A data frame containing the computed scores, with each column representing a score.
#' 
#' @keywords internal
#' 
#' @export
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

#' Fast Calculate Score Profile
#'
#' Calculates the multi-scale KL divergence score profile with dyadic time scales (0, 2, 4, ...).
#'
#' @param W matrix; the affinity matrix.
#' @param max_time integer; the maximum diffusion time. Default is 2^15. The actual maximum diffusion time may be shorter if all nodes converge beforehand.
#' @param init_state matrix; the row-wise normalized initial state matrix.
#' @param P_ls list; the list of diffusion operators. Default is NULL.
#' @param score_ls character vector; the list of scores to compute. See \code{\link{get_score_profile}} for details. Default is \code{c("score0")}.
#'
#' @return A data frame containing the computed scores for different diffusion times.
#'
#' @keywords internal
#'
#' @export
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
    P_ls = ConstructDiffusionOperators(W,max_time)
  }
  
  # Calculate multi-scale KL divergence
  score_df = do.call(cbind,lapply(1:length(P_ls),function(i){
    P = P_ls[[i]]
    state = fastMatMult(init_state, as.matrix(P))
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

#' Fast Calculate Score Profile for Large Data
#'
#' Calculates the multi-scale KL divergence score profile for large datasets with dyadic time scales (0, 2, 4, ...).
#'
#' @param W matrix; the affinity matrix.
#' @param max_time integer; the maximum diffusion time. Default is 2^15. The actual maximum diffusion time may be shorter if all nodes converge beforehand.
#' @param init_state matrix; the row-wise normalized initial state matrix.
#' @param score_ls character vector; the list of scores to compute. See \code{\link{get_score_profile}} for details. Default is \code{c("score0")}.
#'
#' @return A data frame containing the computed scores for different diffusion times.
#'
#' @keywords internal
#'
#' @export
fast_calculate_score_profile_largeData <- function(W, max_time = 2^15, init_state,
                                                   score_ls = c("score0")){
  cat("Calculate LMD score profile for large data...\n")
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  final_state = rep(1/nrow(W),nrow(W)) # bi-stochastic
  
  # Calculate transition matrix
  cat("Run doubly stochastic on affinity matrix...\n")
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
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = floor(log(max_time,2)) + 1,
                                   complete = "=",   # Completion bar character
                                   incomplete = "-", # Incomplete bar character
                                   current = ">",    # Current bar character
                                   clear = FALSE,    # If TRUE, clears the bar when finish
                                   width = 100)
  # state_pre = init_state
  break_flag = FALSE
  for(t in 1:floor(log(max_time,2))){ # t here represents the exponent
    # cat("diffusion time:2^",t,"\n")
    # Automatic decide max_time by checking whether the diag of P -> 1/n
    if(!is_sparse_matrix(P) && (max(abs((ncol(P) * Matrix::diag(P)) * ncol(P) - ncol(P))) < 1e-2 * ncol(P))){
      cat("\nmax diffusion time:2^",t-1,"\n")
      break_flag = TRUE
      break
    }
    # Check the sparsity of P, if too dense, transfer it to dense matrix
    if(is_sparse_matrix(P) && (min(Matrix::rowSums(P == 0)/ncol(P)) < 0.9)){
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


#' Fast Calculate Score Profile with Fine Time Scale
#'
#' Calculates the multi-scale KL divergence score profile with fine time scales (0, 1, 2, 3, ...).
#'
#' @param W matrix; the affinity matrix.
#' @param max_time integer; the maximum diffusion time. Default is 100. The actual maximum diffusion time may be shorter if all nodes converge beforehand.
#' @param init_state matrix; the row-wise normalized initial state matrix.
#' @param score_ls character vector; the list of scores to compute. See \code{\link{get_score_profile}} for details. Default is \code{c("score0")}.
#'
#' @return A data frame containing the computed scores for different diffusion times.
#'
#' @keywords internal
#'
#' @export
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
  cat("Run doubly stochastic on affinity matrix...\n")
  P = Doubly_stochastic(W)
  P = as.matrix(P)
  score_df = get_score_profile(state_0 = init_state,
                               state_t = init_state,
                               state_inf = final_state,
                               # P_t = diag(ncol(P)),
                               score_ls = score_ls)
  colnames(score_df) = paste0(colnames(score_df),"_","0")
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
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

#' Obtain LMD Score (LMDS)
#'
#' Calculates LMD Score (LMDS) from a given score profile.
#'
#' @param score_profile data frame; the score profile containing scores computed at different diffusion times. This is typically the output from \code{fast_calculate_score_profile}, \code{fast_calculate_score_profile_largeData}, or \code{fast_calculate_score_profile_fineTimeScale}.
#' @param correction logical; if TRUE, adjusts the score profile by scaling the score profile by entropy and applying delta correction. Default is FALSE.
#'
#' @return A numeric vector containing the LMDS for each row in the score profile (i.e. each gene).
#'
#' @keywords internal
#'
#' @export
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


#' Fast Get LMDS
#'
#' @param W matrix; the affinity matrix.
#' @param max_time integer; the maximum diffusion time. Default is 2^15. The actual maximum diffusion time may be shorter if all nodes converge beforehand.
#' @param init_state matrix; the row-wise normalized expression matrix.
#' @param P_ls list; the list of diffusion operators. Default is NULL.
#' @param correction logical; if TRUE, adjusts the score profile by scaling the score profile by entropy and applying delta correction. Default is FALSE.
#' @param largeData logical; if TRUE, uses \code{\link{fast_calculate_score_profile_largeData}}. Default is TRUE.
#' @param highres logical; if TRUE, uses \code{\link{fast_calculate_score_profile_highres}}. Default is FALSE.
#'
#' @return A list containing:
#' \item{score_profile}{data frame; the computed score profile for different diffusion times.}
#' \item{cumulative_score}{numeric vector; the LMDS for each row in the score profile.}
#'
#'
#' @export
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

#' Calculate Localized Marker Detector Score (LMDS) for each gene
#'
#' Computes the LMD score of each gene given an expression matrix and cell space.
#'
#' @param expression matrix; the gene by cell expression matrix.
#' @param feature_space matrix; the cell by coordinate matrix (e.g., 20 principal components).
#' @param knn integer; the number of nearest neighbors for constructing the graph. Default is 5.
#' @param kernel logical; if TRUE, uses a Gaussian kernel. Otherwise, uses a kNN binarized graph. Default is FALSE.
#' @param max_time integer; the maximum diffusion time. The actual maximum diffusion time may be shorter if all genes converge beforehand. Default is 2^20.
#' @param adjust_bridge logical; if TRUE, connects disconnected components of the graph using Minimum Spanning Trees. Default is TRUE.
#' @param self_loop integer; the weight for self-connections. Default is 1.
#' @param score_correction logical; if TRUE, adjusts the LMD profile by delta correction. Default is FALSE.
#' @param largeData logical; if TRUE, uses functions optimized for large matrix multiplication. Default is TRUE.
#' @param highres logical; if TRUE, uses fine time scales (0, 1, 2, 3, ...). If FALSE, uses dyadic time scales (0, 2, 4, ...). Default is FALSE.
#' @param min_cell integer; removes genes expressing in fewer than this number of cells. Default is 5.
#'
#' @return A list containing:
#' \item{score_profile}{data frame; the computed score profile for different diffusion times.}
#' \item{cumulative_score}{numeric vector; the LMD score for each gene.}
#'
#'
#' @export
LMD <- function(expression, feature_space, knn = 5, 
                kernel = FALSE, max_time = 2^20, adjust_bridge = TRUE, self_loop = 1,
                score_correction = FALSE, largeData = TRUE, highres = FALSE, min_cell = 5, kernel_used = NULL){
  if(any(colnames(expression) != rownames(feature_space))){stop("Cells in expression mtx and feature space don't match.")}
  if(!kernel){
    W = ConstructKnnGraph(knn = knn, feature_space = feature_space, adjust_disconnection = adjust_bridge, self_loop = self_loop)$'graph'
  }else{
    if(kernel_used == "SNN"){
      W = ConstructSNNGraph(knn = knn, feature_space = feature_space, self_loop = self_loop)$'graph'
    }else if(kernel_used == "Gaussian"){
      W = ConstructGaussianGraph(knn = knn, feature_space = feature_space, alpha = 1, coef = 2, epsilon = 1e-3, self_loop = self_loop)$'graph'
    }else{
      stop("Please provide graph type: SNN or Gaussian")
    }
  }
  rho = RowwiseNormalize(expression)
  rho = rho[,colnames(W),drop = FALSE]
  rho = rho[which(apply(rho,1,function(x) sum(x>0) >= min_cell))
            ,,drop = FALSE] # sanity check & remove genes which express at less than 5 cells
  cat(sprintf("Remove %d genes which express in less than %d cells\n",
              nrow(expression) - nrow(rho),min_cell))
  res = fast_get_lmds(W = W, max_time = max_time,
                      init_state = rho, 
                      correction = score_correction, 
                      largeData = largeData,
                      highres = highres)
  return(res)
}

knee_point = function(Vecs, plot.fig = FALSE){
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
  if(plot.fig == TRUE){
    plot(Vecs,pch = 20, xlab = "Gene Index", ylab = "LMD Score")
    points(knee_point,Vecs[knee_point],col='red',pch = 20)
    cat("knee_point: ", knee_point, "\n")
  }
  return(knee_point)
}

#' Show LMD Result
#'
#' Display the results of LMD.
#'
#' @param res.lmd list; the output of the \code{\link{LMD}}, containing the LMD score profile and LMD scores.
#' @param kneeplot logical; if TRUE, a knee plot of the LMD score distribution will be shown.
#'
#' @return A list containing:
#' \item{gene_table}{data frame; a table with genes, their LMD scores, and their ranks.}
#' \item{gene_rank}{named numeric vector; the rank of each gene.}
#' \item{cut_off_gene}{named numeric vector; the top-ranked genes up to the knee point in the LMD score distribution.}
#'
#'
#' @export
show_result_lmd <- function(res.lmd, kneeplot = FALSE){
  # res.lmd: output of LMD (LMD profile)
  score = res.lmd$cumulative_score
  score = sort(score)
  df = data.frame(score = score)
  df$'rank' = 1:nrow(df)
  # print(head(df,n = 10))
  gene_rank = setNames(df$'rank',rownames(df))
  return(list(gene_table = df, gene_rank = gene_rank, cut_off_gene = gene_rank[1:knee_point(score, plot.fig = kneeplot)]))
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
  cat("Run doubly stochastic on affinity matrix...\n")
  P = Doubly_stochastic(W)
  cat("self-weight:",diag(P)[1],"\n")
  
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = floor(log(max_time,2)) + 1,
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
    if(is_sparse_matrix(P) && (min(Matrix::rowSums(P == 0)/ncol(P)) < 0.9)){
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
