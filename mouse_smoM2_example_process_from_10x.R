###### Preprocess all samples
###### Rihao Qu

source('/data/rq25/seurat_wrapper_funs.R')
root_dir <- "/data/rq25/peggy_10x_mouse_smom2_202112/"
sample_data_dir <- "/data/rq25/peggy_10x_mouse_smom2_202112/data/"
analysis_results_dir <- "/data/rq25/peggy_10x_mouse_smom2_202112/seurat_analysis/result/"
prefix <- "SmoM2_data"

#####
##### Load required packages for analysis
#####
library(Seurat)
library(scales)
library(dplyr)
library(viridis)
library(grid)


#####
##### Data import
#####
setwd(root_dir)
sample_names <- list.dirs(sample_data_dir, full.names = F, recursive = F) #sample_names <- c("KO", "WT") ##if you want to specify samples by yourself
sample_names
data_S_list_v0 <- Load10xData(sample_data_dir, sample_names)


data_S_list_v1 <- data_S_list_v0

##### Quality control
#####
setwd(analysis_results_dir)
get_quality_vlnplot(data_S_list_v1, file = get_output_name("quality_vlnplot.pdf", prefix, "figure"))
get_quality_vlnplot(data_S_list_v1, file = get_output_name("quality_vlnplot_log.pdf", prefix, "figure"), log = T)
data_S_list_v2 <- FilterCells(data_S_list_v1, ngene_lth = 500, ngene_hth = 5000, mt_hth = 10)
metric_report <- get_metric_summary(data_S_list_v1, data_S_list_v2, ngene_lth = 500, ntranscript_lth = 2500)
metric_report
write.csv(metric_report, file = get_output_name("metric_summary.csv", prefix), quote = F)


#####
##### Merge data
#####
if (T){
  data_S_list <- lapply(X = data_S_list_v2, FUN = function(x) {
    x <- NormalizeData(x) 
    #x <- RunALRA(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  })
}

data_S <- MergeData(data_S_list)
table(data_S$orig.ident)

DefaultAssay(data_S) <- "RNA"
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)
data_S <- ScaleData(data_S)

#####
##### Dimensionality reduction
#####
data_S <- RunPCA(data_S, npcs = 30, verbose = F)
plot(data_S@reductions$pca@stdev)

data_S <- RunUMAP(data_S, dims = 1:30)
data_S <- RunUMAP(data_S, dims = 1:30, n.components = 3)
DimPlot(data_S, reduction = "umap", group.by = "orig.ident", shuffle = T)

DimPlot(data_S, reduction = "umap", group.by = "orig.ident", shuffle = T)
DimPlot(data_S, reduction = "umap", group.by = "Phase", shuffle = T) 
DimPlot(data_S, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident")
DimPlot(data_S, reduction = "umap", group.by = "Phase", split.by = "orig.ident") 

data_S <- RunTSNE(
  data_S, tsne.method = "FIt-SNE", check_duplicates = FALSE, do.fast = TRUE, seed.use=3, dims = 1:30, perplexity = 100, #300
  fast_tsne_path="/home/jz437/git/FIt-SNE/bin/fast_tsne", ann_not_vptree=FALSE, nthreads=12
)
#data_S <- RunTSNE(
#  data_S, dims = 1:30,dim.embed = 3, perplexity = 100
#)

DimPlot(data_S, reduction = "tsne", group.by = "orig.ident", shuffle = T)
DimPlot(data_S, reduction = "tsne", group.by = "Phase", shuffle = T) 
DimPlot(data_S, reduction = "tsne", group.by = "orig.ident", split.by = "orig.ident")
DimPlot(data_S, reduction = "tsne", group.by = "Phase", split.by = "orig.ident") 

data_S$log_nCount_RNA <- log(data_S$nCount_RNA)
data_S$log_nFeature_RNA <- log(data_S$nFeature_RNA)

FeaturePlotNew(data_S, features = c("nCount_RNA", "nFeature_RNA"), reduction = "tsne") 
FeaturePlotNew(data_S, features = c("log_nCount_RNA", "log_nFeature_RNA"), reduction = "tsne") 


#####
##### Overlay expression of gene markers
#####
DefaultAssay(data_S) <- "RNA"
markers_to_check <- c("Sox2", "Dkk1", "Dkk2", "Lef1", "Cdkn1a", "Bmp3", "Ptch1", "Bmp4", "Ccnb1")
markers_to_check <- c("Lhx2", "Shh", "Fgf20", "Dkk4", "Col1a1", "Krt10", "Krt14", "Krt5", "Gli1")#c("Sox2", "Dkk1", "Dkk2", "Lef1", "Cdkn1a", "Bmp3", "Ptch1", "Bmp4", "Ccnb1")

FeaturePlot(
  data_S, 
  features = markers_to_check, #"Pax7", 
  ncol = 3,
  #cols = c("gray","red"), 
  #pt.size = 0.5, 
  reduction = "umap", 
  #dims = c(dim1,dim2),
  order =T
) 

gene <- markers_to_check[9]
#gene <- "eYFP"
data_S$tmp.ident <- paste0(data_S$orig.ident, "|", gene)
FeaturePlot(
  data_S, 
  features = gene, #"Pax7", 
  #ncol = 3,
  split.by = "tmp.ident",
  #cols = c("gray","red"), 
  #pt.size = 0.25, 
  reduction = "umap", 
  #dims = c(dim1,dim2),
  order =T
) & theme(legend.position = c(1,0.5),
          plot.margin = unit(rep(0.5,4), "cm")) & 
  scale_color_gradientn(colours = c('gray75', 'blue'), limits = c(0, max(data_S[["RNA"]]@data[gene,]))) #&


FeaturePlot(
  data_S, 
  features = c("Hba-a1", "Hbb-bs", "Hbb-y", "Acta2", "Cox4i2", "Ebf1"), #"Pax7", 
  ncol = 3,
  #cols = c("gray","red"), 
  #pt.size = 0.5, 
  reduction = "umap", 
  #dims = c(dim1,dim2),
  order =T
) 
#####
##### Clustering
#####
data_S <- FindNeighbors(data_S, dims = 1:30)
data_S <- FindClusters(data_S, resolution = 0.3)
data_S$cluster <- data_S$RNA_snn_res.0.3 #integrated
DimPlot(data_S, reduction = "umap", label = T, group.by = "cluster")
DimPlot(data_S, reduction = "umap", label = T, group.by = "cluster", split.by = "orig.ident")
DimPlot(data_S, reduction = "tsne", label = T, group.by = "cluster")
DimPlot(data_S, reduction = "tsne", label = T, group.by = "cluster", split.by = "orig.ident")

#####
##### Find cluster-specific markers
#####

#data_S <- data_S_CTL_B_vs_T[which(rownames(data_S_CTL_B_vs_T) %notin% grep("^TR[ABDG]",rownames(data_S_CTL_B_vs_T),value = TRUE)),]
DefaultAssay(data_S) <- "RNA"
Idents(data_S) <- "cluster"
cluster_markers <- FindAllMarkers(
  data_S, only.pos = T, min.diff.pct = 0.09
)
cluster_markers$pct.diff <- cluster_markers$pct.1 - cluster_markers$pct.2
write.xlsx(
  cluster_markers, file = paste0("./data/",prefix, "_cluster_markers.xlsx")#SMC_
)
saveRDS(cluster_markers, paste0("./data/",prefix, "_cluster_markers.rds"))#SMC_

### Heatmap on cluster-specific markers
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 4, wt = avg_log2FC)
#data_S_scaled_for_heatmap <- FindVariableFeatures(data_S, nfeatures = nrow(data_S))
#data_S_scaled_for_heatmap <- ScaleData(data_S_scaled_for_heatmap)
cells <- sample(colnames(data_S), 20000)
DoHeatmap(data_S, features = top10$gene, cells = cells, slot = "scale.data") + theme(axis.text.y = element_text(size = 8)) + scale_fill_viridis()



data_S <- subset(data_S_peggy_smom2_all, cluster %in% c(0:3, 6:7))
###
###Go back to the initial analysis
###








#############
#####
#####DM
#####
data_S <- subset(data_S_peggy_smom2_dermal, cluster %in% c(2,4) & orig.ident == "E13_SmoM2")
data_S <- subset(data_S_peggy_smom2_dermal, orig.ident == "E13_Cont")

DimPlot(data_S)

DefaultAssay(data_S) <- "RNA"
data_S <- NormalizeData(data_S)
data_S <- FindVariableFeatures(data_S)
data_S <- ScaleData(data_S, vars.to.regress = c("S.Score", "G2M.Score"))

#####
##### Dimensionality reduction
#####
#DefaultAssay(data_S) <- "integrated"
DefaultAssay(data_S)
data_S <- RunPCA(data_S, npcs = 30, verbose = F)
plot(data_S@reductions$pca@stdev)

data_S <- RunUMAP(data_S, dims = 1:30)

DimPlot(data_S, reduction = "umap", group.by = "orig.ident")
DimPlot(data_S, reduction = "umap", group.by = "Phase", split.by = "orig.ident") 
DimPlot(data_S, reduction = "umap", group.by = "cluster")

K <- 10
nPC <- 30
dists <- as.matrix(dist(data_S@reductions$pca@cell.embeddings[,1:nPC]))
sigma_list <- c()
for (i in 1:nrow(dists)){
  sigma_list[i] <- sort(dists[i,])[K]
}

library(FNN)
knn_res <- get.knn(data_S@reductions$pca@cell.embeddings[,1:nPC], k=K)
str(knn_res)
affnity_matrix <- matrix(0, nrow=nrow(dists), ncol=ncol(dists))
for (i in 1:nrow(affnity_matrix)){
  affnity_matrix[i,knn_res$nn.index[i,]] <- exp(-knn_res$nn.dist[i,]^2/(sigma_list[i]^2))
}
dim(affnity_matrix)
sum(affnity_matrix[100,]>0)
#plot(density(affnity_matrix))

affnity_matrix_2 <- (affnity_matrix + t(affnity_matrix))/2

normalized_vec <- sqrt(1/apply(affnity_matrix_2, 1, sum))
#normalized_vec2 <- sqrt(1/apply(affnity_matrix_2, 2, sum))

affnity_matrix_3 <- t(affnity_matrix_2 * normalized_vec) * normalized_vec
#normalized_vec <- sqrt(apply(affnity_matrix_3, 1, sum))
#normalized_vec2 <- sqrt(apply(affnity_matrix_3, 2, sum))
#normalized_vec[1:10]
#normalized_vec2[1:10]

E_list <- rARPACK::eigs_sym(affnity_matrix_3, k = 30)
str(E_list)
plot(1:30, E_list$values[1:30]^10)
#sum(E_list$vector[10,]^2)

diffu_emb <- matrix(0, nrow = nrow(E_list$vector), ncol = 30)
for(i in 1:30){
  diffu_emb[,i] <- E_list$vector[,i]*normalized_vec*(E_list$values[i]^10)
}
colnames(diffu_emb) <- paste0("EV_", 1:ncol(diffu_emb))
rownames(diffu_emb) <- colnames(data_S)

#plot(diffu_emb[,2], diffu_emb[,4])

data_S[["diffu50_tr"]] <- CreateDimReducObject(embeddings = diffu_emb, key = "EV_", assay = DefaultAssay(data_S))
#dplot <- DimPlot(data_S, reduction = "diffu50_tr", dims = c(2,3), group.by = "orig.ident", pt.size = 2.5, shuffle = T, cols = hue_pal()(2)[2:1]) #order = sample(colnames(data_S)), 
#dplot[[1]]$layers[[1]]$aes_params$alpha = .25
#dplot
#data_S <- subset(data_S, orig.ident != "E14M2WT_MYT")

DimPlot(data_S, reduction = "diffu50_tr", dims = c(2,4), group.by = "Phase")#, split.by = "orig.ident")
DimPlot(data_S, reduction = "diffu50_tr", dims = c(2,3), group.by = "Phase", split.by = "orig.ident")
DimPlot(data_S, reduction = "diffu50_tr", dims = c(2,3), group.by = "cluster", split.by = "orig.ident")

DefaultAssay(data_S) <- "RNA"
dim1 <- 2
dim2 <- 3
key <- "diffu50_tr"
markers_to_check <- c("Sox2", "Dkk1", "Dkk2", "Lef1", "Cdkn1a", "Bmp3", "Ptch1", "Bmp4", "Ccnb1")
FeaturePlot(
  data_S, 
  features = markers_to_check, #"Pax7", c( "Sfrp2", "Scx", "Meox1", "Acta2"),#
  ncol = 3,
  #nrow = 3,
  #cols = c("gray","red"), 
  #pt.size = 0.5, 
  reduction = key, 
  dims = c(dim1,dim2),
  #slot = "data", 
  order =T
) & xlim(min(data_S[[key]]@cell.embeddings[,dim1]),max(data_S[[key]]@cell.embeddings[,dim1])) & 
  ylim(min(data_S[[key]]@cell.embeddings[,dim2]),max(data_S[[key]]@cell.embeddings[,dim2]))

