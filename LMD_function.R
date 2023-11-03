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
    stop("Disconnected Components, Please increase knn.")
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
  W = Symmetric_KNN_graph(knn = knn, feature_space = feature_space)$'graph'
  rho = Rowwise_normalize(expression)
  res = fast_calculate_multi_score(W = W, init_state = rho, max_time = max_time)
  return(res)
}

# Visualization =========
FeaturePlot_custom <- function(value,coord,value_name = NULL,title_name = NULL,order_point = TRUE){
  if(length(value)!=nrow(coord)){stop("Unmatched Dimension!")}
  df = cbind(coord[,1:2],value = value)
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
  if(ncol(init_state)!=graph_node_number | nrow(init_state)!=1){stop("Unmatched dimension")}
  check_step = sort(unique(floor(log(check_time,base = 2))))
  check_step = check_step[check_step >= 1]
  include_init = FALSE
  if(0 %in% check_time){
    include_init = TRUE
  }
  check_time = 2^check_step
  if(is.null(P_ls)){P_ls = Obtain_Pls(W,max(check_time))}
  multi_state = sapply(P_ls[check_step],function(P){
    state = fastMatMult(init_state, P)
  })
  if(include_init){
    multi_state = cbind(t(init_state),multi_state)
    check_time = c(0,check_time)
  }
  colnames(multi_state) = check_time
  degree_node = rowSums(W) # Normalize multi_state with degree of node for visualization
  pl = lapply(seq(ncol(multi_state)),function(i){
    legend_break_label = seq(0, max(multi_state[,i]), length.out = 2)
    if(max(legend_break_label)>1e-3){
      legend_break_label = signif(legend_break_label, digits = 2)
    }else{
      legend_break_label = formatC(legend_break_label, format = "e", digits = 0)
    }
    p = suppressWarnings({
      FeaturePlot_custom(multi_state[,i]/max(multi_state[,i]),coord,
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
      ggtitle(plot.title1) + labs(color = "Expression") +
      theme(legend.title = element_text(size = 8),
            panel.grid = element_blank(),
            panel.background = element_blank()) 
    p1
  })
  names(pl) = levels(feature_partition)
  return(pl)
}
Visualize_score_pattern <- function(res.lmd, genes = NULL, label_class = NULL, facet_class = NULL, add_point = NULL){
  score_df = res.lmd$diffusion_KL_score
  score_df["score0_0"] = 0
  if(!all(genes %in% rownames(score_df))){stop("Genes not found!")}
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
  x_breaks = unique(as.numeric(df$step))
  names(x_breaks) = pmax(log(x_breaks,base = 2),0)
  names(x_breaks) = ifelse(names(x_breaks) == "0",names(x_breaks),paste0("2^",names(x_breaks)))
  x_breaks = setNames(names(x_breaks),x_breaks)
  df$step = as.factor(as.numeric(df$step))
  levels(df$step) = x_breaks[levels(df$step)]
  
  if(!is.null(label_class)){
    p = ggplot(data = df, mapping = aes(x = step, y = score, color = label)) + 
      geom_line(aes(group = gene)) + labs(x = "Time", y = "Normalized Diffusion KL Score")  + 
      geom_text(data = df %>% group_by(gene) %>% slice_head(n = 4) %>% slice_tail(n = 1),
                aes(label = gene), hjust = 0) + theme(axis.text.x = element_text(angle = 45,hjust = 1))
  }else{
    p = ggplot(data = df, mapping = aes(x = step, y = score, color = gene)) + 
      geom_line(aes(group = gene)) + labs(x = "Time", y = "Normalized Diffusion KL Score")
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
    p = p + facet_wrap(~facet,scales = "free")
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



