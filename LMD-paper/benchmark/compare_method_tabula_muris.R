cran_packages <- c("readxl","dplyr","BiocManager","RColorBrewer")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)

# Download/Load Data ===========
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")

dir.path <- dir.path0
data_source = "tabular_muris"
folder.path <- file.path(dir.path,data_source)
dir.create(folder.path, recursive=T)
tissue_download_link = c(
  "marrow_facs" = "https://figshare.com/ndownloader/files/13092380",
  "pancreas_facs" = "https://figshare.com/ndownloader/files/13092386",
  "lung_facs" = "https://figshare.com/ndownloader/files/13092194",
  "marrow_droplet" = "https://figshare.com/ndownloader/files/ 13089821"
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
data_source = "tabular_muris"
tissue_name = names(tissue_download_link)[1]

dir.path <- dir.path0
folder.path <- file.path(dir.path,data_source)
tiss <- readRDS(file.path(folder.path,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
n_dim = 20
feature_space = Embeddings(tiss[["pca"]])[,1:n_dim]
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]

# Run Each Method and save results ========
dir.path <- dir.path0
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"run_methods_function.R"))
folder.path = file.path(dir.path,"benchmark",data_source)
dir.create(folder.path, recursive = T)
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot","semi","marcopolo")
# Create RunTime Table
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
  #   RunMarcopolo(tiss, selected_genes, dir.file)
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
dir.path <- dir.path0
folder.path <- file.path(dir.path,data_source,"ground_truth_geneset")
dir.create(folder.path, recursive=T)
file_name = paste0(tissue_name,"_ground_truth_c1.txt")
if(!file.exists(file.path(folder.path,file_name))){
  avg_exp = AverageExpression(subset(tiss,features = selected_genes), assays = "RNA", slot = "counts", group.by = "cell_ontology_class") %>% as.data.frame()
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

## Criterion2: CellMarkerDB =============
file_name = paste0(tissue_name,"_ground_truth_c2.txt")
if(!file.exists(file.path(folder.path,file_name))){
  if(!file.exists(file.path(folder.path,"Cell_marker_All.xlsx"))){
    options(timeout=6000)
    download.file("http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_All.xlsx",
                  destfile = file.path(folder.path,"Cell_marker_All.xlsx"), method = 'libcurl')
  }
  # Load CellMarkerDB
  cell_marker_db = readxl::read_xlsx(file.path(folder.path,"Cell_marker_All.xlsx")) 
  tissue_class = ifelse(grepl("marrow",tissue_name),"Bone marrow",
                        ifelse(grepl("lung",tissue_name),"Lung",
                               ifelse(grepl("pancreas",tissue_name),"Pancreas",NA)))
  cell_marker_db = cell_marker_db %>% filter(species == "Mouse") %>% 
    filter(tissue_class == tissue_class) %>% 
    filter(!is.na(cellontology_id)) %>% filter(!is.na(PMID)) %>% 
    mutate(gene_name1 = marker) %>% mutate(gene_name2 = Symbol) %>% 
    select(c("gene_name1","gene_name2")) 
  # Match gene name
  celldb_marker = union(cell_marker_db$gene_name1,cell_marker_db$gene_name2)
  celldb_marker = unlist( lapply(celldb_marker,function(i){
    id = grep(paste0("^", i, "$"), selected_genes, ignore.case = TRUE)
    if(length(id)==0){return(NA)}
    else{return(selected_genes[id])}
  }) )
  celldb_marker = celldb_marker[!is.na(celldb_marker)] %>% unique()
  write.table(celldb_marker,file = file.path(folder.path, file_name))
}

# Load Rank Table ==========
dir.path <- dir.path0
folder.path.rank <- file.path(dir.path,"benchmark",data_source)
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot","semi","marcopolo")
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
folder.path.gt <- file.path(dir.path,data_source,"ground_truth_geneset")
auc_df = do.call(rbind,lapply(method_ls,function(method){
  vec = quick_marker_benchmark(setNames(df_benchmark[,method],df_benchmark$gene),
                               folder_path = folder.path.gt, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Method = method)
}))
write.table(auc_df,file = file.path(folder.path.rank, paste0(tissue_name,"_auroc.csv")),row.names = FALSE)


# Concordance table ==========
celldb_marker = read.table(file.path(folder.path.gt, paste0(tissue_name,"_ground_truth_c2.txt")))[,1]
df_benchmark = read.table(file.path(folder.path.rank, paste0(tissue_name,"_benchmark_rank_table.csv")),
                          header = TRUE)
# Among top N genes, how many of them are ground-truth?
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot","semi","marcopolo")
unlist(lapply(method_ls,function(method){
  vec = setNames(df_benchmark[,method],df_benchmark$gene)
  n = round(nrow(df_benchmark) * 0.05)
  sum(names(sort(vec))[1:n] %in% celldb_marker)
}))

# Density Index ========
dir.path = dir.path0
for(tissue_name in names(tissue_download_link)[1:3]){
tiss <- readRDS(file.path(dir.path,data_source,paste0(tissue_name,".rds")))
folder.path.rank <- file.path(dir.path,"benchmark",data_source)
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot","semi","marcopolo")
df_benchmark = read.table(file.path(folder.path.rank, paste0(tissue_name,"_benchmark_rank_table.csv")),
                          header = TRUE)
rownames(df_benchmark) = df_benchmark$gene
n_dim = 20
df_density_index =do.call(rbind,lapply(c(method_ls,"All genes"),function(method){
  if(method == "All genes"){
    top_gs = rownames(df_benchmark)
    topn = length(top_gs)
  }else{
    top_gs = rownames(df_benchmark)[order(df_benchmark[,method])]
    topn = c(seq(50,90,10),seq(100,900,100),seq(1000,3000,500),seq(4000,9000,1000),seq(10000,nrow(df_benchmark),2000))
    # topn = c(seq(50,90,10),seq(100,1000,100))
  }
  result = do.call(rbind, lapply(topn, function(top){
    features = top_gs[1:top]
    # remove low-quality cells: cells express in less than 5 genes
    tmp_dat = subset(tiss,features = features)[["RNA"]]@counts
    pass_QC_cell = names(which(colSums(tmp_dat > 0) >= 5))
    subtiss = subset(tiss, cells = pass_QC_cell)
    subtiss <- subtiss %>% ScaleData(features = features) %>% RunPCA(features = features, npcs = n_dim)
    
    # Low-quality_rate
    rate = 1 - length(pass_QC_cell) / ncol(tmp_dat)
    # Density-index
    nn.dist = FNN::get.knn(Embeddings(subtiss, reduction = "pca"), k = 10, algorithm = "kd_tree")$nn.dist
    sdVec <- stats::na.omit(subtiss@reductions$pca@stdev)
    length_scale <- sqrt(sum(sdVec^2)) # length-scale * sqrt(2*(N-1)/N) = d_rms
    s2 = length_scale/mean(nn.dist)
    c(rate, s2)
  }))
  colnames(result) = c("LowQualityRate","DensityIndex")
  result = data.frame(TopGene = topn,Method = method,result)
  result
}) )

write.table(df_density_index,file = file.path(folder.path.rank, paste0(tissue_name,"_DI.csv")),row.names = FALSE)

# error bar & null distribution
topn = c(seq(50,90,10),seq(100,900,100),seq(1000,3000,500),seq(4000,9000,1000),seq(10000,nrow(df_benchmark),2000))
df_density_index_null = do.call(rbind,lapply(topn,function(top){
  set.seed(top)
  seed_ls = sample(1:1e5,size = 100)
  random_top_gs = lapply(seed_ls,function(sed){
    set.seed(sed)
    sample(rownames(df_benchmark),top,replace = FALSE)
  })
  result = unlist(lapply(random_top_gs, function(g_ls){
    features = g_ls
    # remove low-quality cells: cells express in less than 5 genes
    tmp_dat = subset(tiss,features = features)[["RNA"]]@counts
    pass_QC_cell = names(which(colSums(tmp_dat > 0) >= 5))
    subtiss = subset(tiss, cells = pass_QC_cell)
    
    subtiss <- subtiss %>% ScaleData(features = features) %>% RunPCA(features = features, npcs = n_dim)
    # Density-index
    nn.dist = FNN::get.knn(Embeddings(subtiss, reduction = "pca"), k = 10, algorithm = "kd_tree")$nn.dist
    sdVec <- stats::na.omit(subtiss@reductions$pca@stdev)
    length_scale <- sqrt(sum(sdVec^2))
    s2 = length_scale/mean(nn.dist)
    s2
  }))
  result = data.frame(TopGene = top,DensityIndex = result)
}) )
write.table(df_density_index_null,file = file.path(folder.path.rank, paste0(tissue_name,"_DI_null.csv")),row.names = FALSE)
}