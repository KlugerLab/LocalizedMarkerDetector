# =======================
# Visualization =========
# =======================

#' Visualize Graph
#'
#' Visualizes a graph using a given affinity matrix and node layout.
#'
#' @param affinity_m matrix; symmetric affinity matrix of the graph.
#' @param label vector; node labels (optional). Default is NULL.
#' @param layout matrix; a 2D coordinate matrix of cells to visualize the graph.
#'
#' @return A ggplot2 object representing the plotted graph.
#'
#' 
#' 
#' @importFrom igraph graph_from_adjacency_matrix V
#' @import ggraph
#' @import ggplot2
#' @import RColorBrewer
#' 
#' @export
VisualizeGraph <- function(affinity_m, label = NULL, layout){
  layout <- as.matrix(layout)
  g <- igraph::graph_from_adjacency_matrix(affinity_m, 'undirected', weighted = TRUE)
  igraph::V(g)$frame.color <- NA
  
  if(!is.null(label)){
    igraph::V(g)$color <- label
    p = ggraph::ggraph(g,layout = layout) + 
      geom_edge_link(ggplot2::aes(width = weight),color = "grey") + 
      geom_node_point(ggplot2::aes(color = color), size = 1) + 
      scale_edge_width_continuous(range = c(0.5, 1), breaks = c(0.5, 1)) + 
      theme_graph()
    if(class(label) != "numeric"){
      coldef <- setNames(
        colorRampPalette(brewer.pal(12, "Paired"))(length(unique(label))),
        unique(label))
      p = p + scale_color_manual(values = coldef)
    }
  }else{
    p = ggraph::ggraph(g,layout = layout) + 
      geom_edge_link(ggplot2::aes(width = weight),color = "grey") + 
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

#' Custom Feature Plot
#'
#' Creates a custom feature plot for visualizing a feature's values on a reduced dimension embedding (e.g., UMAP or t-SNE).
#'
#' @param coord A Seurat object or a matrix of coordinates. If a Seurat object is provided, embeddings will be extracted from Seurat.
#' @param value A vector of values to plot. If `coord` is a Seurat object, `value` can be a column name from `coord@meta.data`.
#' @param reduction Character; the name of the reduction embedding to use (e.g., "umap", "tsne"). Default is NULL.
#' @param dims A vector specifying which dimensions to use. Default is c(1,2).
#' @param value_name Character; the name of the value to be used in the plot legend. Default is NULL.
#' @param title_name Character; the title of the plot. Default is NULL.
#' @param order_point Logical; whether to order points by value. Default is TRUE.
#'
#' @return A ggplot2 object representing the custom feature plot.
#'
#' 
#' 
#' @import ggplot2
#' @import RColorBrewer
#' 
#' @export
CustomFeaturePlot <- function(coord, value, reduction = NULL, dims = 1:2, value_name = NULL,title_name = NULL,order_point = TRUE){
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
  df[,1] = as.numeric(df[,1]); df[,2] = as.numeric(df[,2]);
  if(order_point){df = df[order(df[,3]),]}
  p <- ggplot(df, aes(x=!!as.name(colnames(coord)[1]), y=!!as.name(colnames(coord)[2]), color=value)) + 
    geom_point(size=0.5)
  if(is.numeric(value)){
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

#' Visualize Diffusion
#'
#' Creates a series of plots visualizing the diffusion process of a gene's expression over time on a given graph.
#'
#' @param coord matrix; 2D coordinates of cells for visualization.
#' @param init_state matrix; row-wise normalized expression matrix for the initial state.
#' @param P_ls list; list of diffusion operators.
#' @param W matrix; symmetric affinity matrix of the graph.
#' @param check_time vector; time points selected for visualization, e.g., \code{c(2, 16, 128)}.
#' @param gene_name character; the name of the gene for visualization.
#' @param gene_color character; the color of expression level for visualization. Default is "blue".
#'
#' @return A ggplot2 object representing a series of plots visualizing the diffusion process.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom patchwork wrap_plots plot_annotation
#' 
#' @export
VisualizeDiffusion <- function(coord, init_state, P_ls = NULL, W = NULL, check_time = NULL, gene_name = NULL, gene_color = "blue"){
  init_state = as.matrix(init_state)
  tmp <- if (!is.null(P_ls) && length(P_ls) > 0) {
    as.matrix(P_ls[[1]])
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
  if(is.null(P_ls)){P_ls = ConstructDiffusionOperators(W,max(check_time))}
  check_time = check_time[check_time %in% names(P_ls)]
  sub = match(check_time,names(P_ls))
  multi_state = sapply(P_ls[sub],function(P){
    P = as.matrix(P)
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
      CustomFeaturePlot(coord = coord, value = multi_state[,i]/max(multi_state[,i]),
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


#' Plot Feature Modules
#'
#' Creates feature plots for visualizing the average expression of gene modules across cells.
#'
#' @param dat A Seurat object or a matrix of gene expression data.
#' @param coord matrix; 2D coordinates of cells for visualization. Required if \code{dat} is not a Seurat object.
#' @param feature_partition factor; a vector indicating the partition of features (genes) into modules.
#' @param reduction character; the name of the reduction embedding to use (e.g., "umap", "tsne"). Default is NULL.
#' @param assays character; the assay to use from the Seurat object. Default is "RNA".
#'
#' @return A list of ggplot2 objects representing the feature plots for each module.
#'
#' 
#' 
#' @import ggplot2
#' @import RColorBrewer
#' 
#' @export
CustomModulePlot <- function(dat, coord = NULL, feature_partition, reduction = NULL, assays = "RNA"){
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
  feature_partition = feature_partition[!is.na(feature_partition)]
  feature_partition = as.factor(feature_partition)
  pl <- lapply(levels(feature_partition), function(level){
    genes = names(feature_partition)[feature_partition == level]
    if(length(genes)==0) {return(NULL)}
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
  pl = Filter(Negate(is.null), pl)
  return(pl)
}

#' Visualize Score Pattern
#'
#' Visualizes the LMD score profiles for selected genes over time.
#'
#' @param score_profile matrix; the LMD score profile, output of \code{fast_get_lmds}.
#' @param genes vector; selected genes for visualization.
#' @param label_class vector; gene class labels for visualizing genes in different colors. Default is NULL.
#' @param facet_class vector; gene facet labels for visualizing genes in different panels. Default is NULL.
#' @param add_point vector; time points to highlight. Default is NULL.
#' @param dyadic logical; if TRUE, represents the axis text in a 2^ format. Default is TRUE.
#' @param text logical; if TRUE, adds gene labels to each curve. Default is FALSE.
#' @param normalize logical; if TRUE, normalizes the LMD profile by the plateau. Default is FALSE.
#'
#' @return A ggplot2 object representing the score patterns for the selected genes.
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import cowplot
#' @import dplyr
#' 
#' @export
VisualizeScorePattern <- function(score_profile, genes = NULL, 
                                    label_class = NULL, facet_class = NULL, 
                                    add_point = NULL, dyadic = TRUE, text = FALSE, normalize = FALSE){
  score_df = score_profile
  if(!all(genes %in% rownames(score_df))){stop("Genes not found!")}
  profiles = names(which(table(sub("_\\d+", "", colnames(score_df)))>1))
  score_df = score_df[genes,]
  if(("score0" %in% profiles) & normalize){
    score_df = score_df/score_df$'max_score0'
  }
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
      geom_line(aes(group = interaction(gene,profiles))) + 
      labs(x = "Time", y = "Normalized Diffusion KL Score") + theme_half_open() + background_grid()
    if(text == TRUE){
      p = p +
        geom_text(data = df %>% group_by(interaction(gene,profiles)) %>% slice_head(n = 4) %>% slice_tail(n = 1),
                  aes(label = gene), hjust = 0) + theme(axis.text.x = element_text(angle = 45,hjust = 1))
    }
    
  }else{
    p = ggplot(data = df, mapping = aes(x = step, y = score, color = gene, linetype = profiles)) + 
      geom_line(aes(group = interaction(gene,profiles))) + 
      labs(x = "Time", y = "Normalized Diffusion KL Score") + theme_half_open() + background_grid()
  }
  if(length(unique(df$profiles)) == 1) {p = p + guides(linetype = "none")}
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

#' Visualize Gene Distance Heatmap
#'
#' Creates a heatmap for visualizing a pairwise distance matrix of genes. The sidebar and side tree represent the partition of genes.
#'
#' @param dist matrix; a pairwise distance matrix of genes.
#' @param gene_partition factor; a vector indicating the partition of genes into modules. Default is NULL.
#' @param gene_hree hclust; a hierarchical clustering object of genes. Default is NULL.
#' @param filtered logical; if TRUE, filters out some genes for each partition based on the SML method (Parisi et al., 2014). Default is TRUE.
#'
#' @return A ggplot2 object representing the heatmap of the gene pairwise distance matrix with gene partitions.
#'
#' @references
#' Parisi, F., Strino, F., Nadler, B., & Kluger, Y. (2014). Ranking and combining multiple predictors without labeled data. Proceedings of the National Academy of Sciences, 111(4), 1253-1258. doi:10.1073/pnas.1219097111
#' 
#' @import ggplot2
#' @import RColorBrewer
#' @import pheatmap
#' 
#' @export
VisualizeGeneHeatmap <- function(dist, gene_partition = NULL, gene_hree = NULL, filtered = TRUE){
  if(is.null(gene_partition)){gene_partition = ClusterGenes(dist, filtered = filtered)}
  if(is.null(gene_hree)){gene_hree = ClusterGenes(dist, return_tree = TRUE, filtered = filtered)[[2]]}
  module_color = setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(nlevels(gene_partition)),
    levels(gene_partition))
  module_color[is.na(names(module_color))] = "lightgrey"
  p = pheatmap(as.matrix(dist),
               annotation_col = data.frame(Module = gene_partition),
               cluster_rows = gene_hree,
               cluster_cols = gene_hree,
               annotation_colors = list(Module = module_color),
               treeheight_row = 0,
               col=colorRampPalette(c("firebrick3", "white"))(99),
               show_colnames = FALSE, 
               show_rownames = FALSE,
               annotation_legend = FALSE,
               # main = paste(distance_method,
               #              clustering_method,sep = "-"),
               silent = TRUE)
  return(ggplotify::as.ggplot(p))
}

#' Visualize Gene t-SNE
#'
#' Creates a t-SNE plot for visualizing gene partitions based on a given distance matrix.
#'
#' @param dist matrix; a pairwise distance matrix of genes.
#' @param gene_partition factor; a vector indicating the partition of genes into modules. Default is NULL.
#' @param filtered logical; if TRUE, filters out some genes for each partition based on the SML method (Parisi et al., 2014). Default is TRUE.
#'
#' @return A ggplot2 object representing the t-SNE visualization of gene partitions.
#'
#' 
#' @references
#' Parisi, F., Strino, F., Nadler, B., & Kluger, Y. (2014). Ranking and combining multiple predictors without labeled data. Proceedings of the National Academy of Sciences, 111(4), 1253-1258. doi:10.1073/pnas.1219097111
#'
#' @import ggplot2
#' @import RColorBrewer
#'
#' @export
VisualizeGeneTSNE <- function(dist, gene_partition = NULL, filtered = TRUE){
  if(is.null(gene_partition)){gene_partition = ClusterGenes(dist, filtered = filtered)}
  set.seed(233)
  df = data.frame(Rtsne::Rtsne(dist, is_distance = TRUE, 
                               dim = 2, perplexity = 30, max.iter = 500)$Y)
  colnames(df) = c("tSNE_1","tSNE_2")
  rownames(df) = labels(dist)
  gene_partition = as.factor(gene_partition)
  
  df[names(gene_partition),"group"] = gene_partition
  module_color = setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(nlevels(gene_partition)),
    levels(gene_partition))
  module_color[is.na(names(module_color))] = "lightgrey"
  centroids <- df %>%
    group_by(group) %>%
    summarise(tSNE_1 = mean(tSNE_1), tSNE_2 = mean(tSNE_2))
  p = ggplot(df, aes(x=tSNE_1, y=tSNE_2, color=group)) + 
    geom_point(size=1) + scale_color_manual(values = module_color) + 
    geom_text(data = centroids, aes(label = group), vjust = -1, color = "black") +
    theme_minimal() + ggtitle("TSNE Visualization of Genes") + labs(color = "Modules")
  return(p)
}

