cran_packages <- c("readxl","dplyr","BiocManager")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)
library(LocalizedMarkerDetector)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"run_methods_function.R"))

# mask the LocalizedMarkerDetector function for testing purpose
findDisconnectedComponents <- function(adj_matrix) {
  # Create a graph from the adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  # Find the connected components
  components <- igraph::components(g)
  # Get the number of disconnected components
  num_components <- components$no
  # Get the corresponding nodes for each component
  component_nodes <- split(igraph::V(g), components$membership)
  list(number_of_components = num_components, components = component_nodes)
}
findShortestEdge <- function(component1, component2, data) {
  # Calculate all pairwise distances between nodes in the two components
  distances <- outer(component1, component2, Vectorize(function(x, y) dist(data[c(x,y), ])))
  # Find the minimum distance and the corresponding nodes
  min_distance_idx <- which(distances == min(distances), arr.ind = TRUE)
  return(list(from = component1[min_distance_idx[1]], to = component2[min_distance_idx[2]], distance = min(distances)))
}
ConstructKnnGraph_from_embeddings <- function(knn = 5, feature_space, self_loop = 1){
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
  filtered_components <- lapply(res$components, function(comp) if (length(comp) >= 10) comp else NULL)
  filtered_components <- Filter(Negate(is.null), filtered_components)
  filtered_node_names <- unlist(lapply(filtered_components, function(comp) names(comp)))
  A = A[rownames(A)%in%filtered_node_names,
        colnames(A)%in%filtered_node_names]
  if(dim(A)[1] < nrow(feature_space)){
    cat("Remove ", nrow(feature_space) - dim(A)[1], " singleton cells.\n")
  }
  
  A = as.matrix(A)
  diag(A) = self_loop
  W = pmin(A + t(A),1) # unweighted
  diag(W) = self_loop
  
  component = findDisconnectedComponents(A)$components
  component_partition <- unlist(lapply(names(component),function(i){
    setNames(rep(i,length(component[[i]])),names(component[[i]]))
  }))
  return(list(graph = W, adj_matrix = A,
              component = component,
              component_partition = component_partition))
}
FindEdges <- function(filtered_components, feature_space){
  cat(length(filtered_components), " Disconnected Components, fully connecting.\n")
  
  # Find shortest edges between components
  edgesToAdd <- lapply(1:(length(filtered_components)-1), function(i) {
    lapply((i + 1):length(filtered_components), function(j) {
      findShortestEdge(filtered_components[[i]], filtered_components[[j]], feature_space)
    })
  })
  edgesToAdd = do.call(rbind.data.frame, lapply(unlist(edgesToAdd, recursive = FALSE),
                                                function(x){data.frame(from = names(x$'from'),
                                                                       to = names(x$'to'),
                                                                       distance = x$distance)} )) %>% distinct()
  edges_group = edgesToAdd;
  edges_group[,1:2] = apply(edges_group[,1:2],2,function(y){
    apply(do.call(rbind,lapply(filtered_components,function(x) y %in% names(x))),2,function(x) which(x))
  })
  
  # Filter edges by MST
  g_meta <- igraph::graph_from_data_frame(edges_group, directed = FALSE)
  igraph::E(g_meta)$weight <- edges_group$distance
  
  mst_edges <- igraph::mst(g_meta, algorithm = "prim")
  mst_edges <- as.data.frame(igraph::get.edgelist(mst_edges))
  edgesToAdd_filtered = edgesToAdd[match(do.call(paste, mst_edges),do.call(paste, edges_group[,1:2])),]
  
  return(list(full_edges = edgesToAdd,
              mst_edges = edgesToAdd_filtered))
}
AddEdges_to_graph <- function(edgesToAdd, A, self_loop = 1){
  # edgesToAdd_filtered = edgesToAdd
  # Add edges
  cat("add ", nrow(edgesToAdd), " edges\n")
  for (i in 1:nrow(edgesToAdd)) {
    A[edgesToAdd[i,1],edgesToAdd[i,2]] = 1
  }
  filtered_components = findDisconnectedComponents(A)$components
  if (length(filtered_components) >= 2) {
    stop("Disconnected Components, Please increase knn or connect them.")
  }
  
  A = as.matrix(A)
  diag(A) = self_loop
  W = pmin(A + t(A),1) # unweighted
  diag(W) = self_loop
  
  return(list(graph = W, adj_matrix = A))
}


# Load Data =======
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")

dir.path <- dir.path0
data_source = "tabular_muris"
tissue_name = "marrow_facs"
tiss <- readRDS(file.path(dir.path,data_source,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]

# Show Cell graph with disconnected components ----
knn = 5
# PCA20 graph
feature_space = Embeddings(tiss[["pca"]])[,1:20]
graph_20pc_ori = ConstructKnnGraph_from_embeddings(knn = knn,
                                                   feature_space = feature_space)
res = graph_20pc_ori
p = VisualizeGraph(affinity_m = res$graph, layout = Embeddings(tiss[["tsne"]])[colnames(res$graph),],
               label = res$component_partition[colnames(res$graph)] ) + ggtitle("Original kNN (20PC)")
ggsave(filename = "~/20pc.png", 
       plot = p, width = 6, height = 5)

# TSNE graph
feature_space = Embeddings(tiss[["tsne"]])[,1:2]
graph_tsne_ori = ConstructKnnGraph_from_embeddings(knn = knn,
                                        feature_space = feature_space)
res = graph_tsne_ori
p = VisualizeGraph(affinity_m = res$graph, layout = Embeddings(tiss[["tsne"]])[colnames(res$graph),],
                    label = res$component_partition[colnames(res$graph)]) + ggtitle("Original kNN (TSNE)")
ggsave(filename = "~/tsne.png", 
       plot = p, width = 6, height = 5)

# Adding Edges
edge_ls = FindEdges(res$component, feature_space)

graph_tsne_fully = AddEdges_to_graph(edge_ls$full_edges,A = graph_tsne_ori$adj_matrix)
graph_tsne_mst = AddEdges_to_graph(edge_ls$mst_edges,A = graph_tsne_ori$adj_matrix)

res = graph_tsne_fully
p = VisualizeGraph(affinity_m = res$graph, 
                    layout = Embeddings(tiss[["tsne"]])[colnames(res$graph),],
                    label = graph_tsne_ori$component_partition[colnames(res$graph)]) + ggtitle("kNN + fully connect")
ggsave(filename = "~/tsne_fully.png", 
       plot = p, width = 6, height = 5)

res = graph_tsne_mst
p = VisualizeGraph(affinity_m = res$graph, 
                   layout = Embeddings(tiss[["tsne"]])[colnames(res$graph),],
                   label = graph_tsne_ori$component_partition[colnames(res$graph)]) + ggtitle("kNN + mst connect")
ggsave(filename = "~/tsne_mst.png", 
       plot = p, width = 6, height = 5)

# Test consistency -----------
res = lapply(list(graph_20pc_ori,
                  graph_tsne_ori,
                  graph_tsne_fully,
                  graph_tsne_mst),function(graph){
                    W = graph$graph
                    rho = RowwiseNormalize(dat)
                    rho = rho[, colnames(W), drop = FALSE]
                    rho = rho[which(apply(rho, 1, function(x) sum(x > 0) >= 5)), 
                              , drop = FALSE]
                    res = fast_get_lmds(W = W, max_time = 2^20, init_state = rho)
                    return(res)
                  })

score_df = do.call(cbind.data.frame,lapply(res,function(x) x$'cumulative_score'))
colnames(score_df) = c("PC20","TSNE","TSNE_fully","TSNE_MST")
score_ranked <- data.frame(score_df) %>% mutate(across(everything(), rank))

## Jaccard Index ----------
base = "PC20"
top_genes = c(50,seq(100,1000,100))
jaccard_index_bw_sets = do.call(cbind,lapply(top_genes,function(top_gene)
{
  top_gene_ls = apply(score_ranked,2,function(x) rownames(score_ranked)[order(x)][1:top_gene])
  apply(top_gene_ls,2,function(x){
    length(intersect(x,top_gene_ls[,base]))/length(union(x,top_gene_ls[,base]))
  })
}))
colnames(jaccard_index_bw_sets) = top_genes
df = reshape2::melt(jaccard_index_bw_sets,value.name = "JaccardIndex")
colnames(df)[1:2] = c("graph","Top_Genes")
df$Top_Genes = as.factor(df$Top_Genes)

## AUROC ----------
folder.path.gt <- file.path(dir.path,data_source,"ground_truth_geneset")
auc_df = do.call(rbind,lapply(colnames(score_ranked),function(method){
  vec = quick_marker_benchmark(setNames(score_ranked[,method],rownames(score_ranked)),
                               folder_path = folder.path.gt, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Parameter = method)
}))

## Plot -------
ggplot(data = df, aes(x = graph, y = JaccardIndex, group = Top_Genes, color = Top_Genes)) +
  geom_line() +
  geom_point() + labs(x = "# of PC used", y = "Jaccard Index\n (with 20PC)", color = "Top Genes") + geom_vline(xintercept = "20", color="red", linetype="dashed") 
# + 
#   scale_y_continuous(limits = c(0,1))

auc_df = auc_df %>% filter(gt_set %in% c("Top100","CellMarkerDB"))
auc_df$gt_set = factor(auc_df$gt_set)
levels(auc_df$gt) = ifelse(grepl("CellMarkerDB",levels(auc_df$gt)),levels(auc_df$gt),"MacFC(Top100)")

ggplot(data = auc_df, aes(x = Parameter, y = AUROC, group = gt_set, color = gt_set)) +
  geom_line() +
  geom_point() + labs(x = "# of PC used", y = "AUROC", color = "Gold-Standard Genesets") + 
  facet_wrap(~gt_set, scales = "free")
# + 
#   scale_y_continuous(limits = c(0,1)) 

# BarPlot
p1 = ggplot(data = df %>% filter(graph %in% c("TSNE_MST")), 
       aes(x = Top_Genes, y = JaccardIndex, group = 1)) +
  geom_line() + labs(x = "Top Genes", y = "Jaccard Index\n (with 20PC)") + ylim(0,1) + 
  theme_minimal() + ggtitle("Jaccard Index of LMD results b/w two graphs")

p2 = ggplot(data = auc_df %>% filter(Parameter %in% c("PC20", "TSNE_MST")) %>% mutate(Parameter = as.factor(Parameter)), 
       aes(x = gt_set, y = AUROC, group = Parameter, fill = Parameter)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "Gold-Standard Geneset", fill = "Graphs") +  
  theme_minimal() + ggtitle("AUROC evaluation of LMD results of two graphs")

ggsave(filename = file.path(figure_path,"rebuttal-connectivity.png"), plot = p1 + p2, width = 10, height = 5, dpi = 300)
