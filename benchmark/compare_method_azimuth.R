cran_packages <- c("readxl","dplyr","BiocManager","RColorBrewer")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)

# Download/Load Data ===========
dir.path0 <- "/banach1/ruiqi/local_marker"
dir.path <- file.path(dir.path0,"LMD_data")
consortium = "azimuth"
folder.path <- file.path(dir.path,consortium)
dir.create(folder.path, recursive=T)
tissue_download_link = c(
  "human_kidney" = "https://seurat.nygenome.org/azimuth/demo_datasets/kidney_demo_stewart.rds",
  "human_bone_marrow" = "https://seurat.nygenome.org/azimuth/demo_datasets/bmcite_demo.rds",
  "human_pancreas" = "https://seurat.nygenome.org/azimuth/demo_datasets/enge.rds",
  "mouse_motor_cortex" = "https://seurat.nygenome.org/azimuth/demo_datasets/allen_mop_2020.rds",
  "human_lung_hlca" = "https://seurat.nygenome.org/hlca_ref_files/ts_opt.rds",
  "human_motor_cortex" = "https://seurat.nygenome.org/azimuth/demo_datasets/allen_m1c_2019_ssv4.rds"
)
for(tissue_name in names(tissue_download_link)){
  file_name = paste0(tissue_name,".rds")
  if (!file.exists(file.path(folder.path,file_name))) {
    options(timeout=6000)
    download.file(tissue_download_link[tissue_name],
                  destfile = file.path(folder.path,file_name), method = 'libcurl')
    # update version
    tiss <- readRDS(file.path(folder.path,file_name))
    tiss <- UpdateSeuratObject(tiss)
    saveRDS(tiss,file = file.path(folder.path,file_name))
  }
}

# Prepare input data ===========
consortium = "azimuth"
tissue_name = names(tissue_download_link)[1]

for(tissue_name in names(tissue_download_link)){
dir.path <- file.path(dir.path0,"LMD_data")
folder.path <- file.path(dir.path,consortium)
tiss <- readRDS(file.path(folder.path,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
# Annotate based on Azimuth (download result from shiny)
predictions <- read.delim(file.path(folder.path,"azimuth_annotate_result",paste0(tissue_name,"_azimuth_pred.tsv")), row.names = 1)
tiss <- AddMetaData(
  object = tiss,
  metadata = predictions)
# projected.umap <- readRDS(file.path(folder.path,"azimuth_annotate_result",paste0(tissue_name,"_azimuth_umap.Rds")))
# tiss <- tiss[, Cells(projected.umap)]
# tiss[['umap.proj']] <- projected.umap

# Only keep cells which the predicted ct > 10
annotation_group = colnames(predictions)[1]
meta.data = tiss@meta.data[,annotation_group]
tiss <- subset(tiss, cells = colnames(tiss)[meta.data %in% names(table(meta.data)[table(meta.data) > 10])])
# Filter-out genes detected in less than 5 cells
tiss <- subset(tiss, features = names(which(rowSums(tiss[["RNA"]]@counts > 0) >= 5)))
# Generate PC embeddings
DefaultAssay(tiss) = "RNA"
tiss = tiss %>% NormalizeData() %>% FindVariableFeatures() %>%
  ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
# Prepare Input
n_dim = 20
feature_space = as.matrix(tiss@reductions$pca@cell.embeddings[,1:n_dim])
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]

# # Save for Marcopolo
# sc <- import("scanpy", convert = FALSE)
# scvi <- import("scvi", convert = FALSE)
# tiss_adata <- sc$AnnData(
#   X = t(raw_count)
# )
# anndata::write_h5ad(tiss_adata, file.path(folder_path,gsub(".rds",".h5ad",file_name)))

# Run Each Method and save results ========
dir.path <- file.path(dir.path0,"LocalMarkerDetector")
source(file.path(dir.path,"benchmark","run_methods_function.R"))
folder.path = file.path(dir.path,"benchmark",consortium)
dir.create(folder.path, recursive = T)
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot","semi")
df_runtime = data.frame("nGenes" = nrow(dat),
                        "nCells" = ncol(dat),
                        "Method" = method_ls,row.names = method_ls)
for(method in method_ls){
  dir.file = file.path(folder.path, paste0(method,"_",tissue_name,".csv"))
  start.time <- Sys.time()
  if(method == "lmd"){
    RunLMD(dat, feature_space, dir.file)
  }
  if(method == "hvg"){
    RunHVG(tiss, selected_genes, dir.file)
  }
  if(method == "wilcox_no_filter"){
    RunSeuratv4(tiss, selected_genes, feature_space, dir.file)
  }
  if(method == "haystack"){
    RunHaystack(dat, feature_space, dir.file)
  }
  if(method == "hotspot"){
    RunHotspot(tiss, selected_genes, feature_space, dir.file)
  }
  if(method == "semi"){
    RunSEMITONES(dat, feature_space, dir.file)
  }
  # if(method == "marcopolo"){
  #   RunMarcopolo(tiss, selected_genes,dir.file)
  # }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  df_runtime[method,"Time"] = as.numeric(time.taken, units = "mins")
}
write.table(df_runtime,file = file.path(folder.path, paste0(tissue_name,"_runtime.csv")))

# Prepare ground_truth Marker ===========
## Criterion1: Fold-change ==============
#' For each gene, sorting cell types by mean expression for each gene 
#' and computing fold changes between consecutive types, 
#' take the maximum value among the N-1 fold change values, 
#' given N cell types.
#' 
dir.path <- file.path(dir.path0,"LMD_data")
folder.path <- file.path(dir.path,consortium,"ground_truth_geneset")
dir.create(folder.path, recursive=T)
file_name = paste0(tissue_name,"_ground_truth_c1.txt")
if(!file.exists(file.path(folder.path,file_name))){
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
  write.table(max_logfc,file = file.path(folder.path, file_name))
}

## Criterion2: GT Marker =============
#' download from https://azimuth.hubmapconsortium.org/references/; 
#' marker genes used to annotate the corresponding cell type level of the reference data

file_name = paste0(tissue_name,"_ground_truth_c2.txt")
if(!file.exists(file.path(folder.path,file_name))){
  # download from https://azimuth.hubmapconsortium.org/references/
  sheet_name = readxl::excel_sheets(file.path(folder.path, "Azimuth_Cell_marker.xlsx"))
  cell_marker_db_ls <- lapply(sheet_name, function(X) readxl::read_excel(file.path(folder.path, "Azimuth_Cell_marker.xlsx"), sheet = X))
  cell_marker_db_ls <- lapply(cell_marker_db_ls, as.data.frame)
  names(cell_marker_db_ls) <- sheet_name
  sheet = sheet_name[grep(tissue_name,sheet_name)]
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
  write.table(celldb_marker,file = file.path(folder.path, file_name))
}

# Load Rank Table ==========
dir.path <- dir.path0
folder.path.rank <- file.path(dir.path,"LMD_data","benchmark",consortium)

df_benchmark = lapply(method_ls, function(method){
  marker = read.table(file.path(folder.path.rank, paste0(method,"_",tissue_name,".csv")),header = TRUE)
  if(!"gene" %in% colnames(marker)){
    marker$'gene' = rownames(marker)
  }
  marker %>% select(gene,rank) %>% distinct()
})
names(df_benchmark) = method_ls
df_benchmark = bind_rows(df_benchmark, .id = "method")
df_benchmark <- tidyr::pivot_wider(df_benchmark, names_from = method, values_from = rank) %>% as.data.frame()
df_benchmark[is.na(df_benchmark)] = nrow(df_benchmark)
write.table(df_benchmark,file = file.path(folder.path.rank, paste0(tissue_name,"_benchmark_rank_table.csv")),row.names = FALSE)

# AUROC ==========
df_benchmark = read.table(file.path(folder.path.rank, paste0(tissue_name,"_benchmark_rank_table.csv")),
                          header = TRUE)
folder.path.gt <- file.path(dir.path,"LMD_data",consortium,"ground_truth_geneset")
auc_df = do.call(rbind,lapply(method_ls,function(method){
  vec = quick_marker_benchmark(setNames(df_benchmark[,method],df_benchmark$gene),
                               folder_path = folder.path.gt, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Method = method)
}))
write.table(auc_df,file = file.path(folder.path.rank, paste0(tissue_name,"_auroc.csv")),row.names = FALSE)

}
