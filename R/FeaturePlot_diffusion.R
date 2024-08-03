#' Visualize diffusion 
#'
#' @param coord matrix; 2D coordinate of cells for visualization
#' @param init_state matrix; Row-wise normalized expression matrix
#' @param P_ls matrix; Diffused transition matrix
#' @param W matrix; Symmetry affinity matrix of kNN graph
#' @param check_time vector; Select several time point for visualization, e.g. c(2,16,128)
#' @param gene_name string; Select one gene for visualization
#' @param gene_color string; The color of expression level for visualization
#' 
#' 
#'
#' @return
#' A ggplot object
#'
#' @export
#'
#' @examples
#'
#'



FeaturePlot_diffusion <- function(coord, init_state, P_ls = NULL, W = NULL, check_time = NULL, gene_name = NULL, gene_color = "blue"){
  # coord: 2D coordinate of cells for visualization
  # init_state: Row-wise normalized expression matrix
  # P_ls: Diffused transition matrix
  # W: Symmetry affinity matrix of kNN graph
  # check_time: Select several time point for visualization, e.g. c(2,16,128)
  # gene_name: Select one gene for visualization
  # gene_color: The color of expression level for visualization
  
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
  if(is.null(P_ls)){P_ls = Obtain_Pls(W,max(check_time))}
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
