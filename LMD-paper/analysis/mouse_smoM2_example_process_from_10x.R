###### Preprocess all samples
# by Rihao Qu
# original file: /data/rq25/peggy_10x_mouse_smom2_202112/seurat_analysis/seurat.R
# original raw data file: /data/rq25/peggy_10x_mouse_smom2_202112/

#####
##### Load required packages & functions for analysis
#####
library(Seurat)
library(scales)
library(dplyr)
library(viridis)
library(grid)
Load10xData <- function(data_dir, sample_names, sub_rna_dir = "filtered_feature_bc_matrix"){
  data_S_list <- list()
  data_dir0 <- data_dir
  if (TRUE){
    for(i in sample_names){
      data_dir <- file.path(data_dir0, i, sub_rna_dir)
      print(data_dir)
      data_S_list[[i]] <- CreateSeuratObject(Read10X(data.dir = data_dir), project = i)
    }
  }
  if (TRUE){
    if (length(grep("^mt",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^mt"
    if (length(grep("^MT",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^MT"
    if (length(grep("^hg19-MT",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^hg19-MT"
    print(pattern_mt)
    for(i in names(data_S_list)){
      data_S_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
        object = data_S_list[[i]], pattern = pattern_mt
      )
    }
  }
  
  data_S_list
}
MergeData <- function(data_S_list){
  sample_names <- names(data_S_list)
  n_sample <- length(sample_names)
  data_S <- merge(
    x = data_S_list[[1]], data_S_list[sample_names[2:n_sample]]
  )
  data_S
}
FilterCells <- function(data_S_list, ngene_lth = NULL, ngene_hth = NULL, mt_hth, ngene_lqth = NULL, ngene_hqth = NULL){
  if (is.null(ngene_lth)){
    for(i in names(data_S_list)){
      data_S_list[[i]]$ngene_lth <- quantile(data_S_list[[i]]$nFeature_RNA,ngene_lqth)
      data_S_list[[i]]$ngene_hth <- quantile(data_S_list[[i]]$nFeature_RNA,ngene_hqth)
      data_S_list[[i]]$mt_hth <- mt_hth
      data_S_list[[i]] <- subset(
        x = data_S_list[[i]], 
        subset = nFeature_RNA > ngene_lth & nFeature_RNA < ngene_hth & percent.mt <= mt_hth
      )
      data_S_list[[i]] <- RenameCells(data_S_list[[i]], add.cell.id = i)
    }
  }
  if (is.null(ngene_lqth)){
    for(i in names(data_S_list)){
      data_S_list[[i]]$ngene_lth <- ngene_lth
      data_S_list[[i]]$ngene_hth <- ngene_hth
      data_S_list[[i]]$mt_hth <- mt_hth
      data_S_list[[i]] <- subset(
        x = data_S_list[[i]], 
        subset = nFeature_RNA > ngene_lth & nFeature_RNA <= ngene_hth & percent.mt <= mt_hth
      )
      data_S_list[[i]] <- RenameCells(data_S_list[[i]], add.cell.id = i)
    }
  }
  
  data_S_list
}

#####
##### Data import
#####
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")
dir.path <- dir.path0
dir.path <- file.path(dir.path,"smom2")
folder.path <- file.path(dir.path,"raw_data","data_10x_output")

sample_names <- list.dirs(folder.path, full.names = F, recursive = F) #sample_names <- c("KO", "WT") ##if you want to specify samples by yourself
sample_names
data_S_list_v0 <- Load10xData(folder.path, sample_names)

#####
##### Quality control
#####
data_S_list_v2 <- FilterCells(data_S_list_v0, ngene_lth = 500, ngene_hth = 5000, mt_hth = 10)


#####
##### Cell Cycle Annotation
#####
if (TRUE){
  data_S_list <- lapply(X = data_S_list_v2, FUN = function(x) {
    x <- x %>% NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    # s.genes <- cc.genes$s.genes
    # g2m.genes <- cc.genes$g2m.genes
    s.genes = gprofiler2::gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    g2m.genes = gprofiler2::gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
    
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  })
}

#####
##### Merge Data
#####
data_S <- MergeData(data_S_list)
table(data_S$orig.ident)
DefaultAssay(data_S) <- "RNA"
data_S <- data_S %>% NormalizeData() %>% 
  FindVariableFeatures() %>% ScaleData()

#####
##### Dimensionality reduction
#####
data_S <- RunPCA(data_S, npcs = 30, verbose = F)
# plot(data_S@reductions$pca@stdev)
data_S <- RunUMAP(data_S, dims = 1:30)
DimPlot(data_S, reduction = "umap", group.by = "orig.ident", shuffle = TRUE)

#####
##### Clustering
#####
data_S <- data_S %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.3)
data_S$cluster <- data_S$RNA_snn_res.0.3 #integrated
DimPlot(data_S, reduction = "umap", label = T, group.by = "cluster")
DimPlot(data_S, reduction = "umap", label = T, group.by = "cluster", split.by = "orig.ident")


#####
##### Extract Dermal Cells
#####
FeaturePlot(data_S,label = TRUE,features = c("Dkk1","Dkk2","Lef1","Ptch1","Sox2","Sox18"),order = T,ncol = 3) & NoAxes()
data_S_dermal <- subset(data_S, cluster %in% c(0:3, 6:7))
data_S_dermal <- data_S_dermal %>% NormalizeData() %>% 
  FindVariableFeatures() %>% ScaleData() %>%
  RunPCA(npcs = 30, verbose = F) %>% 
  RunUMAP(dims = 1:30) %>% FindNeighbors(dims = 1:30) %>% 
  FindClusters(resolution = 0.3)
data_S_dermal$cluster <- data_S_dermal$RNA_snn_res.0.3 #integrated
DimPlot(data_S_dermal, reduction = "umap", label = T, group.by = "cluster")
DimPlot(data_S_dermal, reduction = "umap", label = T, group.by = "orig.ident")

#####
##### Separate into different conditions & Save Results
#####
folder.path <- file.path(dir.path,"raw_data","preprocess_data")
dir.create(folder.path, recursive=T)
FeaturePlot(data_S_dermal, label = TRUE, features = c("eYFP"), order = T) & NoAxes()
sample_ls = unlist(lapply(c("E13.5", "E14.5"),function(i) paste0(i,"_", c("MUT", "CTL"))))
sample_ls = paste0("smom2_dermal_",sample_ls)
assign(paste0("data_S_",sample_ls[1]), subset(data_S_dermal, cluster %in% c(2,4) & orig.ident == "E13_SmoM2"))
assign(paste0("data_S_",sample_ls[2]), subset(data_S_dermal, orig.ident == "E13_Cont"))
assign(paste0("data_S_",sample_ls[3]), subset(data_S_dermal, cluster %in% c(2,4) & orig.ident == "E14_SmoM2"))
assign(paste0("data_S_",sample_ls[4]), subset(data_S_dermal, orig.ident == "E14_Cont"))
lapply(sample_ls,function(sample_name){
  saveRDS(get(paste0("data_S_",sample_name)), 
          file = file.path(folder.path,sprintf("data_S_%s_v0.rds",sample_name)))
})
