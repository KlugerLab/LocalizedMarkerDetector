# ========================
# Mtx operation===========
# ========================

# # Enable C++11
# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
# # Create the C++ function
# # requireNamespace("Rcpp")
# Rcpp::cppFunction('arma::mat fastMatMult(arma::mat A, arma::mat B) {
#   arma::mat C = A * B;
#   return C;
# }', depends="RcppArmadillo")


# Rcpp::sourceCpp(code='
# #include <Rcpp.h>
# using namespace Rcpp;
# NumericVector calcKL(NumericVector p, NumericVector q) {
#   double sum = 0;
#   for (int j = 0; j < p.size(); ++j) {
#     if(p[j] == 0) continue;
#     double val = p[j] * log(p[j] / q[j]);
#     if (NumericVector::is_na(val) || std::isinf(val)) val = 0.0;
#     sum += val;
#   }
#   return NumericVector::create(sum);
# }
# 
# // [[Rcpp::export]]
# NumericVector fastKLMatrix(NumericMatrix x, NumericMatrix y) {
#   int n = x.nrow();
#   NumericVector out(n);
# 
#   for (int i = 0; i < n; ++i) {
#     out[i] = calcKL(x.row(i), y.row(i))[0];
#   }
#   return out;
# }
# 
# // [[Rcpp::export]]
# NumericVector fastKLVector(NumericMatrix x, NumericVector y) {
#   int n = x.nrow();
#   NumericVector out(n);
# 
#   for (int i = 0; i < n; ++i) {
#     out[i] = calcKL(x.row(i), y)[0];
#   }
#   return out;
# }
# ')

#' Fast Matrix Multiplication
#'
#' This function multiplies two matrices using Armadillo.
#'
#' @param A A numeric matrix.
#' @param B A numeric matrix.
#' @return A numeric matrix which is the product of A and B.
#' @export
#' @useDynLib LocalizedMarkerDetector, .registration = TRUE
#' @import Rcpp
#' @importFrom Rcpp sourceCpp
#' @name fastMatMult
#' @keywords internal
#' @examples
#' A <- matrix(1:4, 2, 2)
#' B <- matrix(5:8, 2, 2)
#' fastMatMult(A, B)
NULL

#' Fast KL Divergence Calculation for Matrices
#'
#' This function calculates the KL divergence between rows of two matrices.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#' @return A numeric vector with the KL divergence values.
#' @keywords internal
#' @export
#' @name fastKLMatrix
#' @examples
#' x <- matrix(runif(10), 5, 2)
#' y <- matrix(runif(10), 5, 2)
#' fastKLMatrix(x, y)
NULL

#' Fast KL Divergence Calculation for Matrix and Vector
#'
#' This function calculates the KL divergence between rows of a matrix and a vector.
#'
#' @param x A numeric matrix.
#' @param y A numeric vector.
#' @return A numeric vector with the KL divergence values.
#' @keywords internal
#' @export
#' @name fastKLVector
#' @examples
#' x <- matrix(runif(10), 5, 2)
#' y <- runif(2)
#' fastKLVector(x, y)
NULL


#' Row-wise Normalize a Matrix
#'
#' Normalizes each row of a matrix so that the sum of the values in each row equals 1. Rows with a sum of zero are excluded.
#'
#' @param x matrix; the matrix to be row-wise normalized.
#'
#' @return A numeric matrix with each row normalized so that the sum of the values in each row equals 1.
#'
#'
#' @export
RowwiseNormalize <- function(x){
  # x: matrix to be rowwise normalized
  
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
sinkhorn_knopp_largedata = function(A, niter = 100, tol = 1e-8, verb = FALSE) {
  # code refer to (https://rdrr.io/github/aaamini/nett/src/R/sinkhorn.R)
  delta = Inf
  for (irep in 1:niter) {
    scale_fac <- 1 / sqrt(Matrix::rowSums(A))
    A <- Matrix::.sparseDiagonal(x = scale_fac) %*% A %*% Matrix::.sparseDiagonal(x = scale_fac)
    
    delta = pmax(max(abs(1 - Matrix::rowSums(A))), max(abs(1 - Matrix::colSums(A))))
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
    return(sinkhorn_knopp_largedata(A = W))
  }
}
