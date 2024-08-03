#' Calculate gene pairwise distance 
#'
#' @param dat matrix; expression matrix
#' @param method string; metric for measuring gene-gene distance (option: pearson, euclidean, KL, jaccard, spearman)
#' 
#'
#' @return
#' A distance matrix
#'
#' @export
#'
#' @examples
#'
#'

Calculate_distance <- function(dat, method){
  # dat: expression matrix
  # method: metric for measuring gene-gene distance (option: pearson, euclidean, KL, jaccard, spearman)
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
