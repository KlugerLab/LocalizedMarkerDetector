# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html#Predict_doublets
#' How to define how many PCs based on Elbow plot
# https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
#' The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.
#' The point where the percent change in variation between the consecutive PCs is less than 0.1%.

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(SeuratWrappers)
library(RColorBrewer)
require(scales)
source("/data/ruiqi/local_marker/LocalMarkerDetector/LMD_function.R", echo = F)
setwd("/banach1/ruiqi/Peggy_scdata/tmp")

FindPC = function(srat){
  stdv <- srat[["pca"]]@stdev
  sum.stdv <- sum(srat[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  return(min.pc)
}

# Load Data --------------------
# /data/rq25/peggy_10x_mouse_smom2_202112/seurat_analysis/result/data/
sample_ls = apply(expand.grid(c("E13.5", "E14.5"), c("MUT", "CTL"), stringsAsFactors = FALSE),1,function(x) paste0(x,collapse = "_"))
sample_ls = paste0("smom2_dermal_",sample_ls[c(1,3,2,4)])
for(i in sample_ls){
  assign(paste0("data_S_",i), 
         readRDS(sprintf("/data/rq25/peggy_10x_mouse_smom2_202112/seurat_analysis/result/data/data_S_peggy_smom2_dermal_%s.rds",i)))
}
# Preprocess ------------
for(i in length(sample_ls)){
  data_S <- get(paste0("data_S_",sample_ls[i]))
  data_S <- data_S %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
  min.pc = FindPC(srat = data_S)
  data_S <- RunUMAP(data_S, dims = 1:min.pc, seed.use = 42)
  data_S <- RunALRA(data_S); DefaultAssay(data_S) <- "RNA"
  assign(paste0("data_S_",sample_ls[i]), data_S)
}
rm(data_S)

# RunLMD ---------------
folder_path = "/banach1/ruiqi/Peggy_scdata"
sample_ls = paste0("smom2_dermal_",sample_ls[c(1,3,2,4)])
sample_name = sample_ls[1]
assign(sprintf("data_S_%s",sample_name),
       readRDS(file.path(folder_path,sprintf("/tmp/data_S_%s.rds",sample_name))))
data_S <- get(sprintf("data_S_%s",sample_name))

dat = as.matrix(data_S[["RNA"]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = names(Gene_detected_count)[(Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)]
dat = dat[selected_genes,,drop = FALSE]
feature_space = as.matrix(data_S@reductions$pca@cell.embeddings[,1:FindPC(srat = data_S)])
tic()
lmd_result = LMD(dat, feature_space, max_time = 2^20, highres = FALSE)
toc()
local_gene = show_result_lmd(lmd_result)

# Group_Localized_genes ----------------
# data_S <- RunALRA(data_S, assay = "RNA")
dat_alra = as.matrix(data_S[["alra"]]@data)[names(local_gene$cut_off_gene),]
distance_method = "jaccard"
dist = Calculate_distance(dat_alra, method = distance_method)
clustering_method = "average"
gene_hree = hclust(dist, method = clustering_method)
gene_partition = dynamicTreeCut::cutreeDynamic(
  dendro = gene_hree, 
  distM = as.matrix(dist), deepSplit = 2,
  pamStage = TRUE, pamRespectsDendro = TRUE,
  minClusterSize = 10)
names(gene_partition) = labels(dist)
gene_partition = as.factor(gene_partition)
module_order = order(sapply(levels(gene_partition), function(mod) {
  median(which(gene_partition[gene_hree$order] == mod))
}))
levels(gene_partition)[module_order] = seq(nlevels(gene_partition))
gene_partition = factor(gene_partition,levels = seq(nlevels(gene_partition)))
rm(module_order,gene_hree,dat_alra,dist)
local_gene$'gene_partition' = gene_partition

saveRDS(local_gene,file = file.path(folder_path,sprintf("/tmp/local_gene_%s.rds",sample_name)))

pl = lapply(sample_ls, function(sample_name){
  DimPlot(get(paste0("data_S_",sample_name)),group.by = "Phase") + ggtitle(sprintf("%s\n%d (%d) genes, %d cells",gsub("smom2_dermal_","",sample_name),nrow(get(paste0("data_S_",sample_name))),length(get(paste0("local_gene_",sample_name))$'gene_rank'),ncol(get(paste0("data_S_",sample_name)))))
})

# Calculate Module score ----------
## Load Seurat Object & Localized Genes
sample_ls = apply(expand.grid(c("E13.5", "E14.5"), c("MUT", "CTL"), stringsAsFactors = FALSE),1,function(x) paste0(x,collapse = "_"))
sample_ls = paste0("smom2_dermal_",sample_ls[c(1,3,2,4)])
sample_name = sample_ls[1]
assign(sprintf("data_S_%s",sample_name),
       readRDS(file.path(folder_path,sprintf("/tmp/data_S_%s.rds",sample_name))))
assign(sprintf("local_gene_%s",sample_name),
       readRDS(file.path(folder_path,sprintf("/tmp/local_gene_%s.rds",sample_name))))

## Calculate Module Activity Score
local_gene <- get(sprintf("local_gene_%s",sample_name))
data_S <- get(sprintf("data_S_%s",sample_name))
data_S = AddModuleActivityScore(data_S, local_gene$gene_partition, 
                                do_local_smooth = TRUE)

# Visualize Gene/Cell module -------------
module_df = data.frame(module = levels(local_gene$gene_partition), 
                       min_rank = unlist(lapply(levels(local_gene$gene_partition),function(i) min(local_gene$gene_rank[names(local_gene$gene_partition)[local_gene$gene_partition == i]]))),
                       cell_num = colSums(data_S@meta.data[,paste0("Module",levels(local_gene$gene_partition))] > 0.5),
                       gene_num = as.numeric(table(local_gene$gene_partition)) )
## Visualize gene avg expression ---------
pl_module = FeaturePlot_meta(data_S,feature_partition = local_gene$gene_partition)
(wrap_plots(pl_module,ncol = 6) & 
    labs(color = "AvgExp") & 
    NoAxes())

## Visualize module activity score (0~1) ---------
modules = module_df$module[module_df$gene_num >= 10 & 
                             (module_df$min_rank < 100 | 
                                module_df$cell_num > 5)]
pl_module = lapply(modules,function(module){
  FeaturePlot(data_S, features = paste0("Module",module), order = TRUE) + 
    # scale_color_gradientn(colours = rev(brewer_pal(palette = "RdYlBu")(10)), limits = c(0,1)) +
    scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) +
    ggtitle(sprintf("Module %s (%d genes)",module, sum(local_gene$gene_partition == module)))
})
names(pl_module) = modules
(wrap_plots(pl_module,ncol = 6) & 
    labs(color = "ModuleScore") & 
    NoAxes()) + plot_layout(guides = "collect")

## Top1 Genes for each module -----------
pl_module = lapply(modules,function(module){
  g = names(local_gene$cut_off_gene)[local_gene$gene_partition == module][1]
  FeaturePlot(data_S, features = g, order = TRUE) + 
    ggtitle(sprintf("%s (Module %s; rank %d)",g, module, 
                    local_gene$cut_off_gene[g]))
})
wrap_plots(pl_module,ncol = 6) & NoAxes()

## Top10 Genes for module i ----------
i = 1
FeaturePlot(data_S,features = names(local_gene$gene_partition)[
  local_gene$gene_partition == i][1:10],order = TRUE,ncol = 5) & 
  NoLegend() & NoAxes()
message(paste(names(local_gene$gene_partition)[
  local_gene$gene_partition == i], collapse = ", "))
for(i in levels(local_gene$gene_partition)){
  cat("Module",i,"\n")
  message(paste(names(local_gene$gene_partition)[
    local_gene$gene_partition == i], collapse = ", "))
}


# Compare LMD Rank -------------
sample_ls = apply(expand.grid(c("E13.5", "E14.5"), c("MUT", "CTL"), stringsAsFactors = FALSE),1,function(x) paste0(x,collapse = "_"))
df_rank = data.frame(row.names = 
                       Reduce(union,lapply(paste0("local_gene_",sample_ls),
                                           function(x) names(get(x)$gene_rank))))
for(i in sample_ls){
  df_rank[names(get(paste0("local_gene_",i))$gene_rank),paste0(i,"_rank")] = get(paste0("local_gene_",i))$gene_rank
}
df_rank[is.na(df_rank)] = nrow(df_rank)
for(i in sample_ls){
  df_rank[names(get(paste0("local_gene_",i))$gene_partition),paste0(i,"_module")] = get(paste0("local_gene_",i))$gene_partition
}
coldef = setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(modules)),modules)
coldef[setdiff(levels(df_rank$EXB_partition),names(coldef))] = "grey"
# scatterplot
partition = "EXB_module"
ggplot(df_rank, aes(x=EXB_rank, y=WTB_rank, color = get(partition))) + geom_point(size = 1)  + 
  geom_abline(intercept = 0, slope = 1, col = "grey", linetype = "dashed") + 
  scale_x_log10(breaks = 10^(0:4), labels = expression(10^0, 10^1, 10^2, 10^3, 10^4)) +
  scale_y_log10(breaks = 10^(0:4), labels = expression(10^0, 10^1, 10^2, 10^3, 10^4)) + 
  expand_limits(x = c(1, 10^4), y = c(1, 10^4)) + 
  scale_color_manual(values = coldef) + labs(color = partition)
# boxplot
long_df_rank = tidyr::pivot_longer(df_rank, 
                                   cols = 1:length(sample_ls), 
                                   names_to = "Rank_Type", 
                                   values_to = "Rank")
long_df_rank$Rank_Type = factor(long_df_rank$Rank_Type,levels = paste0(c("EXB","EX3","WTB","CTL"),"_rank"))
ggplot(long_df_rank %>% filter(get(partition) %in% modules), aes(x=get(partition), y=Rank, color = Rank_Type)) + 
  geom_boxplot(position = position_dodge(0.75), width = 0.7) + labs(x = partition) + scale_y_continuous(breaks = c(seq(0,1000,100),seq(2000,5000,1000),seq(6000,max(long_df_rank$Rank),5000)),   # Ensure all categories are represented
                                                                                                        labels = c(seq(0,1000,100),seq(2000,5000,1000),seq(6000,max(long_df_rank$Rank),5000))) + 
  geom_point(position = position_jitterdodge(), size = 0.5) +
  scale_y_log10(breaks = 10^(0:4), labels = expression(10^0, 10^1, 10^2, 10^3, 10^4))


# Overlaid several modules on MUT / CTL --------
## Define which modules
sample_id = sample_ls[1]; modules = c(16,18)
## Overlaid modules on different samples
for(i in sample_ls[1:2]){ # i: the samples to be overlaid on
  data_S = get(paste0("data_S_",i))
  gene_partition = get(paste0("local_gene_",sample_id))$gene_partition[get(paste0("local_gene_",sample_id))$gene_partition %in% modules]
  data_S = AddModuleActivityScore(data_S, gene_partition = gene_partition, do_local_smooth = FALSE)
  sub = grep(paste0("Module",modules,collapse = "|"),colnames(data_S@meta.data))
  colnames(data_S@meta.data)[sub] = paste0(sample_id,"_",colnames(data_S@meta.data)[sub])
  assign(paste0("data_S_",i),data_S)
  rm(data_S, gene_partition, sub)
}
## Visualize Module Activity Score
data_S = get(paste0("data_S_",i))
local_gene = get(sprintf("local_gene_%s",sample_id))
FeaturePlot(data_S, features = paste0(sample_id,"_Module",modules), order = TRUE)  & NoAxes() & 
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) & NoLegend()

# Overlaid the union of modules on different samples--------
sample_id = sample_ls[1]; modules = c(2,3,4) #1
sample_id = sample_ls[2]; modules = c(3:4) # c(5,7)
gene_partition = get(paste0("local_gene_",sample_id))$gene_partition[get(paste0("local_gene_",sample_id))$gene_partition %in% modules]
module_name = "EX3-G1"
gene_partition = as.factor(setNames(rep(module_name,length(gene_partition)),names(gene_partition)))
for(i in sample_ls[1:2]){
  data_S_tmp = get(paste0("data_S_",i))
  data_S_tmp = AddModuleActivityScore(data_S_tmp, gene_partition = gene_partition, 
                                      do_local_smooth = FALSE)
  assign(paste0("data_S_",i),data_S_tmp)
  rm(data_S_tmp)
}
data_S = AddModuleActivityScore(data_S, gene_partition = gene_partition, do_local_smooth = FALSE)
pl = lapply(c(paste0(c("EXB","WTB"),"_Dermal"),"all"),function(i){
  if(i != "all"){
    data_S = get(paste0("data_S_",i))
  }
  FeaturePlot(data_S, features = paste0("Module",module_name), order = TRUE) + 
    scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) +
    ggtitle(i)
})
(wrap_plots(pl,nrow = 1) & NoAxes()) + plot_layout(guides = "collect")

# Gene-Gene Correlation of one module -------
sample_id = sample_ls[1]
local_gene = get(sprintf("local_gene_%s",sample_id))
genes = names(local_gene$gene_partition)[local_gene$gene_partition == 19]
data_S = get(sprintf("data_S_%s",sample_id))
dat_alra = as.matrix(data_S[["alra"]]@data)[genes,]
dist = Calculate_distance(dat_alra, method = "jaccard")
as.ggplot(pheatmap(1 - as.matrix(dist), 
           treeheight_col = 0, 
           treeheight_row = 0, 
           cluster_rows=FALSE, cluster_cols=FALSE,
           annotation_colors = ann_colors, 
           col=rev(colorRampPalette(c("red","white","blue4"))(99)), 
           show_colnames = FALSE, 
           show_rownames = FALSE, 
           annotation_legend = TRUE,
           legend_breaks = c(seq(0,1,0.2),1),
           legend_labels = c(seq(0,1,0.2),"jaccard index\n"),
           silent = TRUE)) + theme(plot.margin=unit(c(0.1,0,0,0), 'inches'))

# Cell Cycle embedding -------  
## Obtain CC genes & Re-embed PCA-------
### Obtain CC genes
library(gprofiler2)
s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

# local_gene_tmp1 = get(paste0("local_gene_",sample_ls[1]))
# local_gene_tmp2 = get(paste0("local_gene_",sample_ls[2]))
# genes_CC = union(names(local_gene_tmp1$gene_partition)[local_gene_tmp1$gene_partition %in% c(11:13)],
#                  names(local_gene_tmp2$gene_partition)[local_gene_tmp2$gene_partition %in% c(4:7)])
genes_CC = c(s.genes,g2m.genes)

### Merge MUT & CTL
data_S = merge(get(paste0("data_S_",sample_ls[1])),
               get(paste0("data_S_",sample_ls[2])))
data_S$'type' = as.factor(data_S$old.ident)
levels(data_S$'type') = c("CTL","MUT")

data_S <- data_S %>% FindVariableFeatures()
genes = genes_CC
# genes = VariableFeatures(data_S)
# genes = rownames(data_S)
data_S <- data_S  %>%
  ScaleData(verbose = FALSE, features = genes) %>% RunPCA(npcs = 50, verbose = FALSE, features = genes)
data_S <- RunUMAP(data_S, dims = 1:FindPC(srat = data_S), seed.use = 42)
# data_S <- data_S %>% FindNeighbors() %>% FindClusters(res = 0.3)
DimPlot(data_S, group.by = "Phase", split.by = "type",ncol = 1)
DimPlot(data_S, group.by = "type", shuffle = TRUE) + 
  DimPlot(data_S, group.by = "Phase", shuffle = TRUE) 

## Obtain DM ------
data_S <- GeneTrajectory::RunDM(data_S, 
                                reduction = "pca",
                                dims = 1:FindPC(data_S))
DimPlot(data_S, reduction = "dm",group.by = "Phase", dims = c(1:2))  
gene_visual = data.frame(Embeddings(data_S@reductions$dm))
plotly::plot_ly(gene_visual, 
                x = gene_visual[,1], 
                y = gene_visual[,2], 
                z = gene_visual[,3], size = 0.5)
FeaturePlot(data_S, reduction = "dm", 
            features = paste0(sample_id,"_Module",modules), 
            order = TRUE, split.by = "type") & 
  xlim(-0.004,0.002) & ylim(-0.0035,0.004) & 
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1))

## Infer Pseudo-time based on Angle --------
dm1 = Embeddings(data_S@reductions$dm)[,1]; dm2 = Embeddings(data_S@reductions$dm)[,2]
angles = atan2(dm2 - mean(dm2), dm1 - mean(dm1))
angles[angles < 0] = angles[angles < 0] + 2*pi
data_S$'angles' = angles
FeaturePlot(data_S, reduction = "dm", 
            features = "angles", 
            order = TRUE) & 
  xlim(-0.004,0.002) & ylim(-0.0035,0.004) & 
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10)))
cell_order = names(sort(angles))
cell_order = lapply(c("E13_SmoM2","E13_Cont"),function(i){
  cell_order[grepl(i,cell_order)]
})
names(cell_order) = c("MUT","CTL")
## Smooth gene expression vec based on cell order -------
### Extract genes from target gene modules
sample_id = sample_ls[1]; modules = 19
local_gene = get(sprintf("local_gene_%s",sample_id))
genes = names(local_gene$gene_partition)[local_gene$gene_partition %in% modules]
rm(local_gene)
### Obtain normalized expression vector
rho = lapply(setNames(levels(data_S$type),levels(data_S$type)), function(i){
  Rowwise_normalize(subset(data_S, type == i, features = genes)@assays$alra@data)
})
### Smoothing
rho_smoothed = lapply(setNames(levels(data_S$type),levels(data_S$type)), function(condition){
  mtx = do.call(rbind,lapply(genes,function(g){
    if(!g %in% rownames(rho[[condition]])){return(numeric(dim(rho[[condition]])[2]))}
    vec = rho[[condition]][g,cell_order[[condition]]]
    # loess_fit <- loess(vec ~ seq(1:length(vec)), span = 0.1)
    # smoothed_vec <- predict(loess_fit, data.frame(seq(1:length(vec))))
    # smoothed_vec[smoothed_vec < 0] <- 0
    smoothed_vec <- zoo::rollapply(vec, width = 10, FUN = mean, partial = TRUE, align = 'center')
    # plot(vec, type = 'l', col = 'gray', lty = 2,
    #      xlab = 'Cell Order', ylab = 'Normalized Expression Level', main = g)
    # lines(smoothed_vec, col = 'red', lwd = 2)
    smoothed_vec
  }))
  mtx = Matrix(mtx, sparse = TRUE)
  rownames(mtx) = genes; colnames(mtx) = cell_order[[condition]]
  mtx
})
# ks_dist = lapply(setNames(levels(data_S$type),levels(data_S$type)), function(condition){
#   ks_dist = sapply(1:nrow(rho_smoothed[[condition]]),function(i){
#     cat(i,"\n")
#     r = rho_smoothed[[condition]][i,]
#     vec = sapply(i:nrow(rho_smoothed[[condition]]),function(j){
#       r1 = rho_smoothed[[condition]][j,]
#       ks.test(x = r, y = r1)$'statistic'
#     })
#     vec = c(rep(NA,i-1), vec)
#     vec
#   })
#   ks_dist[upper.tri(ks_dist)] = t(ks_dist)[upper.tri(ks_dist)]
#   ks_dist = Matrix(ks_dist,sparse = TRUE)
#   rownames(ks_dist) = colnames(ks_dist) = rownames(rho_smoothed[[condition]])
#   ks_dist
# })

### Obtain jaccard index & Visualize
jaccard_idx = lapply(setNames(levels(data_S$type),levels(data_S$type)), function(condition){
  1 - as.matrix(Calculate_distance(as.matrix(rho_smoothed[[condition]]), method = "jaccard"))
})

as.ggplot(pheatmap(jaccard_idx[["CTL"]], 
                   treeheight_col = 0, 
                   treeheight_row = 0, 
                   cluster_rows=FALSE, cluster_cols=FALSE,
                   annotation_colors = ann_colors, 
                   col=rev(colorRampPalette(c("red","white","blue4"))(100)), 
                   show_colnames = FALSE, 
                   show_rownames = FALSE, 
                   annotation_legend = TRUE,
                   breaks = seq(0,1,1/99),
                   legend_breaks = c(seq(0,1,0.2),1),
                   legend_labels = c(seq(0,1,0.2),"jaccard index\n"),
                   silent = TRUE)) + theme(plot.margin=unit(c(0.1,0,0,0), 'inches')) + 
  ggtitle("CTL")


# Check DA-cluster ----------
# devtools::install_github("KlugerLab/DAseq")
library(DAseq)
# python2use <- "/path/to/your/python"
labels.1 = levels(data_S$type)[1]
labels.2 = levels(data_S$type)[2]
da_cells <- getDAcells(
  X = Embeddings(data_S@reductions$pca)[,1:FindPC(data_S)],
  cell.labels = as.character(data_S$type),
  labels.1 = labels.1,
  labels.2 = labels.2,
  k.vector = seq(50, 500, 50),
  plot.embedding = Embeddings(data_S@reductions$umap)
)
da_cells$pred.plot
da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(-0.7,0.7),
  plot.embedding = Embeddings(data_S@reductions$umap)
)
da_cells$da.cells.plot
data_S$'da.pred' = da_cells$da.pred
da.cell = rep("NA",length = ncol(data_S)); 
da.cell[da_cells$da.up] = paste0(labels.1,"_abundant")
da.cell[da_cells$da.down] = paste0(labels.2,"_abundant")
data_S$'da.cell' = factor(da.cell,levels = c(paste0(c(labels.1,labels.2),"_abundant"),"NA"))
DimPlot(data_S, group.by = "da.cell", shuffle = TRUE, cols = c(hue_pal()(2),"lightgrey")) + 
  ggtitle("DA cells (|DA_score| > 0.7)")
FeaturePlot(data_S, features = "da.pred", ncol = 1) + 
  ggtitle("DA score") & 
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10)))

## Find DA cluster -------
data_S$'DA' = data_S$da.cell
levels(data_S$DA)[3] = "Other"
data_S$'cluster' = 1
data_S <- get_DA_cluster(data_S,
                         resolution = 0.1,
                         partition.by = "cluster")  
DimPlot(data_S, group.by = "DA_cluster", order = T, 
        label = T, label.size = 5, cols = hue_pal()(length(unique(data_S$DA_cluster))), na.value = "lightgray")
Idents(data_S) = data_S$DA_cluster