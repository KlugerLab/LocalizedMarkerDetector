cran_packages <- c("igraph","ggplot2", "cowplot", "RColorBrewer", 
                   "data.table", "dplyr", "patchwork", 
                   "pheatmap", "ggplotify", "ggraph", 
                   "ClusterR", "Rcpp", "RcppArmadillo", "tictoc")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
lapply(cran_packages, require, character.only = TRUE)

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

# Build Cell Graph & Transition Mtx =======
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
Symmetric_KNN_graph <- function(knn = 5, feature_space, adjust = FALSE){
  knn_list <- FNN::get.knn(feature_space, k = knn, algorithm = "kd_tree")
  # Adjacency matrix
  A <- matrix(0,nrow = nrow(feature_space), ncol = nrow(feature_space), dimnames = list(rownames(feature_space), rownames(feature_space) ))
  for(i in 1:nrow(A)){
    A[i, knn_list[[1]][i,] ] = 1
  }
  
  res = findDisconnectedComponents(A)
  filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 5) comp else NULL)
  filtered_components <- Filter(Negate(is.null), filtered_components)
  
  if(adjust){
    while(length(filtered_components) >= 2){
      knn = knn + 1
      knn_list <- FNN::get.knn(feature_space, k = knn, algorithm = "kd_tree")
      # Adjacency matrix
      A <- matrix(0,nrow = nrow(feature_space), ncol = nrow(feature_space), dimnames = list(rownames(feature_space), rownames(feature_space) ))
      for(i in 1:nrow(A)){
        A[i, knn_list[[1]][i,] ] = 1
      }
      
      res = findDisconnectedComponents(A)
      filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 5) comp else NULL)
      filtered_components <- Filter(Negate(is.null), filtered_components)
    }
    print(paste0("Final kNN: ",knn))
  }
  if (length(filtered_components) >= 2) {
    warning("Disconnected Components, Please increase knn.")
  }
  
  filtered_node_names <- unlist(lapply(filtered_components, function(comp) names(comp)))
  A = A[rownames(A)%in%filtered_node_names,
        colnames(A)%in%filtered_node_names]
  
  # Graph affinity matrix
  W = (A + t(A))/2
  diag(A) = 1
  return(list(graph = W,adj_matrix = A, component = filtered_components))
}
Obtain_Pls <- function(W, max_time){
  P = Doubly_stochastic(W)
  # P = Rowwise_normalize(W)
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
plot_knn_graph <- function(affinity_m, label, layout){
  layout <- as.matrix(layout)
  g <- graph_from_adjacency_matrix(affinity_m, 'undirected')
  V(g)$frame.color <- NA
  
  # coldef <- c(brewer.pal(6,"Blues"),brewer.pal(6,"Reds"))
  
  coldef <- setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(length(unique(label))),
    unique(label))
  
  col <- coldef[label]
  names(col) <- label
  
  plot.igraph(g, layout = layout, vertex.size=1,
              vertex.color=col,vertex.label=NA)
  legend("topright",names(coldef),col=coldef,
         pch=16,cex=0.5,bty='n')
}

# Calculating Diffusion score =========
fast_calculate_multi_score <- function(W, max_time = 2^10, init_state, P_ls = NULL){
  if((ncol(W) != ncol(init_state)) & (ncol(W) == nrow(init_state))){
    init_state = t(init_state)
  }
  if(ncol(W) != ncol(init_state)){
    stop("Check the dimension!")
  }
  # degree = rowSums(W)/sum(rowSums(W))
  degree = rep(1/nrow(W),nrow(W))
  
  # Calculate transition matrix
  if(is.null(P_ls)){
    P_ls = Obtain_Pls(W,max_time)
  }
  
  # Calculate multi-scale KL divergence
  score_df = do.call(cbind,lapply(P_ls,function(P){
    state = fastMatMult(init_state, P)
    c(fastKLMatrix(init_state, state),fastKLVector(state,degree))
  }) )
  score_df = cbind(score_df[1:nrow(init_state),],
                   score_df[(nrow(init_state)+1):nrow(score_df),])
  colnames(score_df) = paste0(rep(c("score0_","score1_"), each = ncol(score_df)/2),
                              colnames(score_df))
  score_df = cbind(score_df,fastKLVector(init_state,degree))
  colnames(score_df)[ncol(score_df)] = "score1_0"
  score_df = score_df[,c(1:floor(ncol(score_df)/2),ncol(score_df),(floor(ncol(score_df)/2)+1):(ncol(score_df)-1))]
  rownames(score_df) = rownames(init_state)
  
  sub_score0 = grep("score0",colnames(score_df))
  sub_score1 = grep("score1_0",colnames(score_df))
  score_df = data.frame(score_df[,sub_score0]/score_df[,sub_score1])
  cumulative_score = rowSums(score_df)
  
  return(list(diffusion_KL_score = score_df,cumulative_score = cumulative_score))
}
show_result_lmd <- function(res.lmd, n = 10){
  score = res.lmd$cumulative_score
  score = sort(score)
  df = data.frame(score = score)
  head(df,n = n)
}
LMD <- function(expression, feature_space, knn = 5, max_time = 2^10){
  if(any(colnames(expression) != rownames(feature_space))){stop("Cells in expression mtx and feature space don't match.")}
  W = Symmetric_KNN_graph(knn = 5, feature_space = feature_space)[[1]]
  rho = Rowwise_normalize(expression)
  res = fast_calculate_multi_score(W = W, init_state = rho)
  return(res)
}

# Visualization =========
FeaturePlot_custom <- function(value,coord,value_name = NULL,title_name = NULL){
  if(length(value)!=nrow(coord)){stop("Unmatched Dimension!")}
  df = cbind(coord[,1:2],value = value)
  df = df[order(df[,3]),]
  p <- ggplot(df, aes(x=!!as.name(colnames(coord)[1]), y=!!as.name(colnames(coord)[2]), color=value)) + 
    geom_point(size=0.2)
  if(!is.factor(df[,3])){
    p <- p + scale_color_gradient(low = "lightgrey", high = "blue")
  }
  p = p + theme(legend.title = element_text(size = 8),
          panel.grid = element_blank(),
          panel.background = element_blank()) + labs(color = value_name) + ggtitle(title_name)
  return(p)
}
FeaturePlot_diffusion <- function(coord, init_state, P_ls = NULL, W = NULL, check_time = NULL, gene_name = NULL){
  init_state = as.matrix(init_state)
  tmp <- if (!is.null(P_ls) && length(P_ls) > 0) {
    P_ls[[1]]
  } else if (!is.null(W)) {
    W
  } else {
    NULL
  }
  if(is.null(tmp)){stop("Cell Graph 'W' not found")}
  graph_node_number = ncol(tmp); rm(tmp)
  if((ncol(init_state)!=graph_node_number) & (nrow(init_state) == graph_node_number)){
    init_state = t(init_state)
  }
  if(ncol(init_state)!=graph_node_number | nrow(init_state)!=1){stop("Unmatched dimension")}
  check_step = log(check_time,base = 2)
  if(is.null(P_ls)){P_ls = Obtain_Pls(W,max(check_time))}
  multi_state = sapply(P_ls[check_step],function(P){
    state = fastMatMult(init_state, P)
  })
  colnames(multi_state) = check_time
  pl = lapply(seq(ncol(multi_state)),function(i){
    FeaturePlot_custom(multi_state[,i],coord,value_name = "Density",title_name = paste0("T = ",colnames(multi_state)[i]))
  })
  p = plot_grid(plotlist = pl,ncol = length(check_time)) + labs(title = gene_name)
  title = ggdraw() + draw_label(gene_name, fontface='bold')
  p = plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1))
  return(p)
}
FeaturePlot_meta <- function(dat, coord, feature_partition){
  feature_partition = as.factor(feature_partition)
  pl <- lapply(levels(feature_partition), function(level){
    genes = names(feature_partition)[feature_partition == level]
    plot.title1 <- sprintf("Module %s (%d genes)",level,length(genes))
    df = cbind(coord[colnames(dat),,drop = FALSE],
               value = apply(dat[genes,,drop = FALSE],2,mean))
    p1 <- ggplot(df[order(df$value),], aes(x=!!as.name(colnames(df)[1]), y=!!as.name(colnames(df)[2]), color=value)) + 
      geom_point(size = 0.2) + 
      scale_color_gradient(low = "lightgrey", high = "blue") + 
      ggtitle(plot.title1) +  
      theme(legend.title = element_text(size = 8),
            panel.grid = element_blank(),
            panel.background = element_blank()) 
    p1
  })
  names(pl) = levels(feature_partition)
  return(pl)
}
Visualize_score_pattern <- function(res.lmd, genes = NULL, label_class = NULL, facet_class = NULL){
  score_df = res.lmd$diffusion_KL_score
  score_df["score0_0"] = 0
  if(!all(genes %in% rownames(score_df))){stop("Genes not found!")}
  score_df = score_df[genes,]
  score_df$'gene' = rownames(score_df)
  score_df$'label' = NA
  score_df$'facet' = NA
  
  if(!is.null(label_class)){
    if(length(label_class)!=length(genes)){
      warning("Labels doesn't match genes, ignore labels")
    }else{
      score_df$'label' = label_class
    }
  }
  
  if(!is.null(facet_class)){
    if(length(facet_class)!=length(genes)){
      warning("Facet labels doesn't match genes, ignore facet labels")
    }else{
      score_df$'facet' = facet_class
    }
  }

  df <- reshape2::melt(score_df,
                       id = colnames(score_df)[!grepl("score",colnames(score_df))],
                       variable.name = c("step"), value.name = "score")
  df$step = gsub("score0_","",df$step)
  df$step = as.factor(as.numeric(df$step))
  
  p = ggplot(data = df, mapping = aes(x = step, y = score, color = label)) + 
    geom_line(aes(group = gene)) + labs(x = "Time", y = "Scaled Diffusion KL Score") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)
          ) + 
    geom_text(data = df %>% group_by(gene) %>% slice_head(n = 4) %>% slice_tail(n = 1),
              aes(label = gene), hjust = 0) + theme_minimal_grid()
  if(!is.null(facet_class)){
    p = p + facet_wrap(~facet,scales = "free")
  }
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

## Calculate Cell Module-activity score
Obtain_cell_partition <- function(expr_dat, gene_partition, 
                                  cell_kNN_graph, major_vote = 5){
  # Fit GMM
  if(sum(!names(gene_partition) %in% rownames(expr_dat))){
    stop("gene name doesn't match data")}
  
  meta_g = sapply(levels(gene_partition),function(topic){
    colMeans(expr_dat[names(gene_partition)[gene_partition == topic],,drop = FALSE])
  })
  rownames(meta_g) = colnames(expr_dat)
  
  tic()
  cell_block = apply(meta_g,2,function(vec){GMM_partition(vec)})
  toc()
  # cell_block = as.data.frame(meta_g) %>% reframe(across(everything(), GMM_partition))
  cell_block = as.matrix(cell_block) + 0
  
  # Local smooth
  cell_block = cell_kNN_graph %*% cell_block
  cell_block = (cell_block > major_vote) + 0
  
  colnames(cell_block) = levels(gene_partition)
  rownames(cell_block) = colnames(expr_dat)
  return(cell_block)
}
GMM_partition <- function(vector){
  opt = ClusterR::Optimal_Clusters_GMM(as.matrix(vector),max_clusters = 10, criterion = "BIC", km_iter = 10, em_iter = 5, seed = 1, plot_data = FALSE)
  opt_num = which.min(opt)
  gmm = ClusterR::GMM(as.matrix(vector),gaussian_comps = opt_num, km_iter = 10, em_iter = 5, seed = 1)
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



