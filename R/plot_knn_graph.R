#' Plot a KNN graph
#'
#' @param affinity_m matrix; Symmetry affinity matrix of kNN graph
#' @param label string; Node label
#' @param layout matrix; 2D coordinate of cells for visualization
#'
#' @return
#' A ggplot object
#'
#' @export
#'
#' @examples
#'
#'


plot_knn_graph <- function(affinity_m, label = NULL, layout){
  # affinity_m: Symmetry affinity matrix of kNN graph
  # label: Node label
  # layout: 2D coordinate of cells for visualization
  
  layout <- as.matrix(layout)
  g <- graph_from_adjacency_matrix(affinity_m, 'undirected', weighted = TRUE)
  V(g)$frame.color <- NA
  
  if(!is.null(label)){
    V(g)$color <- label
    p = ggraph(g,layout = layout) + 
      geom_edge_link(aes(width = weight),color = "grey") + 
      geom_node_point(aes(color = color), size = 1) + 
      scale_edge_width_continuous(range = c(0.5, 1), breaks = c(0.5, 1)) + 
      theme_graph()
    if(class(label) != "numeric"){
      coldef <- setNames(
        colorRampPalette(brewer.pal(12, "Paired"))(length(unique(label))),
        unique(label))
      p = p + scale_color_manual(values = coldef)
    }
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