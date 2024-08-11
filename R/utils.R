# =======================================
# Miscellaneous Test Function =============
# =======================================
# Extrinsic Function: FindPC

#' Determine Number of Principal Components
#'
#' Defines the number of principal components (PCs) to retain based on the elbow plot criteria.
#' The function uses two criteria: 
#' 1) The point where the principal components only contribute 5% of the standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
#' 2) The point where the percent change in variation between consecutive PCs is less than 0.1%.
#'
#' @param srat A Seurat object containing the PCA reduction.
#' @param reduction Character; the name of the reduction technique to use. Default is "pca".
#'
#' @return Integer; the number of principal components to retain.
#'
#' @references
#' For more details on the elbow plot method, refer to the [Harvard Chan Bioinformatics Core scRNA-seq training](https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html).
#'
#' @keywords internal
#' @export
FindPC = function(srat, reduction = "pca"){
  stdv <- srat[[reduction]]@stdev
  sum.stdv <- sum(srat[[reduction]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  return(min.pc)
}


generate_pseudo_gene <- function(data){
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
