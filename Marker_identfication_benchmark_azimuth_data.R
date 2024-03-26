library("Seurat")
library("dplyr")
library(RColorBrewer)
library(reticulate)
library(anndata)

setwd("/data/ruiqi/local_marker/dataset/azimuth")
# Load Data ----------
tissue_download_ls = c(
  "human_motor_cortex" = "https://seurat.nygenome.org/azimuth/demo_datasets/allen_m1c_2019_ssv4.rds",
  "human_kidney" = "https://seurat.nygenome.org/azimuth/demo_datasets/kidney_demo_stewart.rds",
  "human_bone_marrow" = "https://seurat.nygenome.org/azimuth/demo_datasets/bmcite_demo.rds",
  "human_pancreas" = "https://seurat.nygenome.org/azimuth/demo_datasets/enge.rds",
  "mouse_motor_cortex" = "https://seurat.nygenome.org/azimuth/demo_datasets/allen_mop_2020.rds",
  "human_lung_hlca" = "https://seurat.nygenome.org/hlca_ref_files/ts_opt.rds"
)
for(tissue_name in names(tissue_download_ls)){
  cat(tissue_name,"\n")
folder_path = "/data/ruiqi/local_marker/dataset/azimuth"
tissue_name = paste0(tissue_name,"_azi_demo")

file_name = paste0(tissue_name,".rds")
options(timeout=6000)
if (!file.exists(file.path(folder_path,file_name))) {
  options(timeout=6000)
  download.file(tissue_download_ls[tissue_name],
                destfile = file.path(folder_path,file_name), method = 'libcurl')
}
tiss <- readRDS(file.path(folder_path,file_name))

# Preprocess ----------------
# Annotate based on Azimuth (download result from shiny)
folder_path = "/data/ruiqi/local_marker/dataset/azimuth/azimuth_annotate_result"
predictions <- read.delim(file.path(folder_path,paste0(tissue_name,"_azimuth_pred.tsv")), row.names = 1)
tiss <- AddMetaData(
  object = tiss,
  metadata = predictions)
projected.umap <- readRDS(file.path(folder_path,paste0(tissue_name,"_azimuth_umap.Rds")))
tiss <- tiss[, Cells(projected.umap)]
tiss[['umap.proj']] <- projected.umap

annotation_group = colnames(predictions)[1]
meta.data = tiss@meta.data[,annotation_group]
tiss <- subset(tiss, cells = colnames(tiss)[meta.data %in% names(table(meta.data)[table(meta.data) > 10])])
DefaultAssay(tiss) = "RNA"
tiss = tiss %>% NormalizeData() %>% FindVariableFeatures() %>%
  ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
n_dim = 30
# tiss <- tiss %>% RunTSNE(dims = 1:n_dim, reduction = "pca")
# DimPlot(tiss, group.by = annotation_group, reduction = "tsne",label = TRUE, cols = coldef)

DefaultAssay(tiss) <- "RNA"
feature_space = as.matrix(tiss@reductions$pca@cell.embeddings[,1:n_dim])
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]
raw_count = as.matrix(tiss@assays$RNA@counts)[selected_genes,,drop = FALSE]
total_count = tiss@meta.data$nCount_RNA

# # Save for Marcopolo
# sc <- import("scanpy", convert = FALSE)
# scvi <- import("scvi", convert = FALSE)
# tiss_adata <- sc$AnnData(
#   X = t(raw_count)
# )
# anndata::write_h5ad(tiss_adata, file.path(folder_path,gsub(".rds",".h5ad",file_name)))

# GroundTruth ----------
folder_path = "./azimuth_benchmark_result"
if (!dir.exists(folder_path)) {
  # If it doesn't exist, create it
  dir.create(folder_path)
}
# Criterion1
if(!file.exists(file.path(folder_path, paste0(tissue_name,"_ground_truth_c1.txt")))){
  avg_exp = AverageExpression(subset(tiss,features = selected_genes), assays = "RNA", slot = "counts", group.by = annotation_group) %>% as.data.frame()
  avg_exp_ordered <- avg_exp %>%
    rowwise() %>%
    mutate(cell_order_val = list(sort(c_across(everything()), decreasing = FALSE))) %>%
    ungroup()
  avg_exp_ordered = avg_exp_ordered %>%
    rowwise() %>%
    mutate(log_fold_change = list(c(NA, diff(log2(unlist(cell_order_val)+1))))) %>%
    ungroup()
  max_logfc = unlist( lapply(avg_exp_ordered$log_fold_change,function(vec){max(vec,na.rm = TRUE)}) )
  names(max_logfc) = rownames(avg_exp)
  max_logfc = sort(max_logfc,decreasing = TRUE)
  write.table(max_logfc,file = file.path(folder_path, paste0(tissue_name,"_ground_truth_c1.txt")))
}

# Criterion2 - GT cellmarker
if(!file.exists(file.path(folder_path, paste0(tissue_name,"_ground_truth_c2.txt")))){
  # download from https://azimuth.hubmapconsortium.org/references/
  sheet_name = readxl::excel_sheets("/data/ruiqi/local_marker/dataset/Cell_Marker/Azimuth_Cell_marker.xlsx")
  cell_marker_db_ls <- lapply(sheet_name, function(X) readxl::read_excel("/data/ruiqi/local_marker/dataset/Cell_Marker/Azimuth_Cell_marker.xlsx", sheet = X))
  cell_marker_db_ls <- lapply(cell_marker_db_ls, as.data.frame)
  names(cell_marker_db_ls) <- sheet_name
  sheet = sheet_name[grep(gsub("_azi_demo","",tissue_name),sheet_name)]
  ct_retrieve = lapply(unique(tiss@meta.data[,annotation_group]),
                       function(ct) {
                         x = cell_marker_db_ls[[sheet]][,"Label"]
                         x[grep(paste0("^", ct, "$"),x,ignore.case = TRUE)]
                       } )
  names(ct_retrieve) = unique(tiss@meta.data[,annotation_group])
  cell_marker_db = cell_marker_db_ls[[sheet]] %>% filter(Label %in% unlist(ct_retrieve))
  celldb_marker = unlist(lapply(cell_marker_db$Markers, function(x) strsplit(x,split = ",")))
  celldb_marker = trimws(celldb_marker)
  celldb_marker = celldb_marker[celldb_marker %in% selected_genes] %>% unique()
  write.table(celldb_marker,file = file.path(folder_path, paste0(tissue_name,"_ground_truth_c2.txt")))
}

# Run Each Method and save results ---------
method_ls = c("lmd","hvg","wilcox_default",
              "haystack","semi","hotspot", "wilcox_no_filter")
source("/data/ruiqi/local_marker/LocalMarkerDetector/LMD_function.R", echo = F)
time_df = data.frame(row.names = method_ls)
for(method in method_ls){
  cat(method,"\n")
  start <- Sys.time()
  if(method == "lmd"){
    res = LMD(dat,feature_space,max_time = 2^20,knn = 5)
    marker = show_result_lmd(res)$'gene_table'
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = TRUE)
  }
  if(method == "hvg"){
    DefaultAssay(tiss) <- "RNA"
    # tiss <- NormalizeData(tiss, normalization.method = "LogNormalize", scale.factor = 1e6)
    marker <- subset(tiss,features = selected_genes) %>% ScaleData(features = selected_genes) %>% 
      FindVariableFeatures(nfeatures = 5000) %>% VariableFeatures() %>% as.data.frame() %>% setNames(c("gene")) %>% mutate(rank = 1:n())
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = FALSE)
  }
  if(method == "haystack"){
    library(singleCellHaystack)
    set.seed(123)
    res <- haystack(x = feature_space, expression = dat)
    marker <- show_result_haystack(res.haystack = res)
    marker[order(marker$log.p.adj),"rank"] = 1:nrow(marker)
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = TRUE)
  }
  
  if(method == "semi"){
    # pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip
    reticulate::py_run_string(
      "
import SEMITONES
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from SEMITONES.cell_selection import from_knn_dist
from SEMITONES.cell_selection import from_gui
from SEMITONES.cell_selection import from_2D_embedding
from SEMITONES.cell_selection import get_cells_from_gui
from sklearn.metrics.pairwise import pairwise_kernels
from SEMITONES.enrichment_scoring import calculate_escores
from SEMITONES.enrichment_scoring import permute
from SEMITONES.enrichment_scoring import sig_interval
from SEMITONES.support_funcs import pairwise_similarities
from SEMITONES.support_funcs import sig_dictionary
from SEMITONES.support_funcs import sig_bool
from SEMITONES.enrichment_scoring import pvals_per_cell

# Load data
data = pd.DataFrame(r.dat.T)  # counts
PC_embed = pd.DataFrame(r.feature_space)
# Reference Cell
g = 8.6e-4  
S = pairwise_kernels(PC_embed, metric='rbf', gamma=g) 
median = np.median(S, axis=0)
start = int(np.argmin(median))
num_refcell = int(data.shape[1]*0.005)
dd_rcells = from_knn_dist(X=PC_embed, 
                          n_ret=num_refcell,
                          start=start, 
                          metric='rbf',
                          metric_params={'gamma': g})
S = pairwise_similarities(PC_embed,
                          query=dd_rcells,
                          metric='rbf',
                          metric_params={'gamma': g})
dd_rcells = data.index[dd_rcells]
escores = calculate_escores(data.values, query=dd_rcells, S=S)
P = permute(data.values)
pscores = calculate_escores(P, query=dd_rcells, S=S)
pvals = pvals_per_cell(escores,pscores,ret = 'q')
result = pd.DataFrame({'escore': escores.abs().max(1), 'padj': pvals.min(1)})
      ")
    
    marker = reticulate::py$result
    rownames(marker) = rownames(dat)
    marker = marker %>% arrange(.,padj,desc(escore))
    marker$'rank' = 1:nrow(marker)
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = TRUE)
  }
  if(method == "hotspot"){
    # pip install git+https://github.com/yoseflab/Hotspot.git
    cell_id = colnames(raw_count)
    gene_id = rownames(raw_count)
    
    reticulate::py_run_string("
import hotspot
import scanpy as sc
import pandas as pd
from scipy.sparse import csc_matrix
adata = sc.AnnData(X = r.raw_count.T,obs = pd.DataFrame(index = r.cell_id), var = pd.DataFrame(index = r.gene_id))
adata.obs['total_counts'] = r.total_count
adata.layers['counts'] = adata.X.copy()
adata.obsm['X_pca'] = r.feature_space
adata.layers['counts_csc'] = csc_matrix(adata.layers['counts'])

hs = hotspot.Hotspot(
    adata,
    layer_key='counts_csc',
    model='danb',
    latent_obsm_key='X_pca',
    umi_counts_obs_key='total_counts'
)
hs.create_knn_graph(weighted_graph=False, n_neighbors=30)

# Determining informative genes
hs_results = hs.compute_autocorrelations()
    ")
    marker = reticulate::py$hs_results
    marker[order(marker$FDR),"rank"] = 1:nrow(marker)
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = TRUE)
  }
  if(method == "wilcox_no_filter"){
    tiss <- FindNeighbors(tiss, dims = 1:n_dim)
    tiss <- FindClusters(tiss, resolution = 1.5, algorithm = 3)
    # marker <- FindAllMarkers(tiss, features = selected_genes, only.pos = FALSE, test.use = "wilcox", logfc.threshold = 0, min.pct = 0)
    marker <- presto:::wilcoxauc.Seurat(X = subset(tiss,features = selected_genes),
                                        assay = "data", 
                                        seurat_assay = "RNA")
    colnames(marker)[1] = "gene"
    gene_rank = marker %>% arrange(.,padj,desc(logFC)) %>% distinct(.,gene, .keep_all = TRUE) %>% select(gene)
    gene_rank$'rank' = 1:nrow(gene_rank)
    marker <- merge(marker,gene_rank,by = "gene")
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = TRUE)
  }
  if(method == "wilcox_default"){
    tiss <- FindNeighbors(tiss, dims = 1:n_dim)
    tiss <- FindClusters(tiss, resolution = 1.5, algorithm = 3)
    marker <- FindAllMarkers(tiss, features = selected_genes, only.pos = FALSE, test.use = "wilcox", logfc.threshold = 0.25,min.pct = 0.1)
    gene_rank = marker %>% arrange(.,p_val_adj,desc(avg_log2FC)) %>% distinct(.,gene, .keep_all = TRUE) %>% select(gene)
    gene_rank$'rank' = 1:nrow(gene_rank)
    marker <- merge(marker,gene_rank,by = "gene")
    write.table(marker,file = file.path(folder_path, paste0(method,"_",tissue_name,".csv")),row.names = TRUE)
  }
  end <- Sys.time()
  
  time_df[method,"Time"] = end - start
}
write.table(time_df, file = file.path(folder_path, paste0(tissue_name,"_runtime.csv")),row.names = TRUE)

# BenchMark AUROC ------------
folder_path = "./azimuth_benchmark_result"
# tissue_name = names(tissue_download_ls)
# tissue_name = paste0(tissue_name,"_azi_demo")
df_benchmark = lapply(method_ls, function(method){
  marker = read.table(file.path(folder_path, paste0(method,"_",tissue_name,".csv")),header = TRUE)
  if(!"gene" %in% colnames(marker)){
    marker$'gene' = rownames(marker)
  }
  marker %>% select(gene,rank) %>% distinct()
})
names(df_benchmark) = method_ls
df_benchmark = bind_rows(df_benchmark, .id = "method")
df_benchmark <- tidyr::pivot_wider(df_benchmark, names_from = method, values_from = rank) %>% as.data.frame()
df_benchmark[is.na(df_benchmark)] = nrow(df_benchmark)
# write.table(df_benchmark,file = file.path(folder_path, paste0(tissue_name,"_benchmark_rank_table.csv")),row.names = FALSE)

auc_df = do.call(rbind,lapply(method_ls,function(method){
  vec = quick_marker_benchmark(setNames(df_benchmark[,method],df_benchmark$gene),
                               folder_path = folder_path, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Method = method)
}))
write.table(auc_df,file = file.path(folder_path, paste0(tissue_name,"_auroc.csv")),row.names = FALSE)

}

tissue_ls = paste0(names(tissue_download_ls),"_azi_demo")
df_auc = do.call(rbind,lapply(tissue_ls,function(tissue_name){df = read.table(file.path(folder_path,paste0(tissue_name,"_auroc.csv")),header = TRUE); df$'Tissue' = tissue_name;df}))
run_time = do.call(rbind,lapply(tissue_ls,function(tissue_name){
  df = read.table(file.path(folder_path,paste0(tissue_name,"_runtime.csv")),header = TRUE); 
  df$'Tissue' = tissue_name;
  df$'Method' = rownames(df);df
}))
