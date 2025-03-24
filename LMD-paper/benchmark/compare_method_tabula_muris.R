cran_packages <- c("readxl","dplyr","BiocManager","RColorBrewer")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"run_methods_function.R"))

# Download/Load Data ===========
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")

dir.path <- dir.path0
data_source = "tabular_muris"
folder.path <- file.path(dir.path,data_source)
if(!dir.exists(folder.path)){
  dir.create(folder.path, recursive=T)
}

tissue_download_link = c(
  "marrow_facs" = "https://figshare.com/ndownloader/files/13092380",
  "pancreas_facs" = "https://figshare.com/ndownloader/files/13092386",
  "lung_facs" = "https://figshare.com/ndownloader/files/13092194",
  "marrow_droplet" = "https://figshare.com/ndownloader/files/13089821"
)
tissue_ls = names(tissue_download_link)
for(tissue_name in tissue_ls){
  file_name = paste0(tissue_name,".rds")
  if (!file.exists(file.path(folder.path,file_name))) {
    options(timeout=6000)
    download.file(tissue_download_link[tissue_name],
                  destfile = file.path(folder.path,file_name), method = 'libcurl')
    # update version
    load(file.path(folder.path,file_name))
    tiss <- UpdateSeuratObject(tiss)
    saveRDS(tiss,file = file.path(folder.path,file_name))
  }
}

# Obtain some subset data
## bone marrow granulocyte
tissue_name = "marrow_facs_granulocyte"
if(!file.exists(file.path(folder.path,paste0(tissue_name,".rds")))){
  tiss <- readRDS(file.path(folder.path,"marrow_facs.rds"))
  tiss <- subset(tiss, cell_ontology_class == "granulocyte")
  
  # Filter-out genes detected in less than 5 cells
  tiss <- subset(tiss, features = names(which(rowSums(tiss[["RNA"]]@counts > 0) >= 5)))
  
  tiss <- tiss %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = 20, verbose = FALSE)
  tiss = RunUMAP(tiss, dims = 1:20)
  saveRDS(tiss,file = file.path(folder.path,paste0(tissue_name,".rds")))
}
tissue_ls = c(tissue_ls, tissue_name)

# Preprocess data ===========
data_S_ls = lapply(tissue_ls, function(tissue_name){
  dir.path <- dir.path0
  folder.path <- file.path(dir.path,data_source)
  tiss <- readRDS(file.path(folder.path,paste0(tissue_name,".rds")))

  # Only keep cells which the predicted ct > 10
  tiss$'annotation_group' = tiss$'cell_ontology_class'
  meta.data = tiss$annotation_group
  # cat(min(table(meta.data)),"\n")
  tiss <- subset(tiss, cells = colnames(tiss)[meta.data %in% names(table(meta.data)[table(meta.data) > 10])])
  # saveRDS(tiss, file.path(folder.path,paste0(tissue_name,".rds")))
  tiss
})
names(data_S_ls) = tissue_ls

# Run ===========
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot",
              "cosg","scmarker","semi","marcopolo")

lapply(tissue_ls, function(tissue_name){
tiss = data_S_ls[[tissue_name]]
# Prepare input data --------
DefaultAssay(tiss) <- "RNA"
n_dim = 20
feature_space = Embeddings(tiss[["pca"]])[,1:n_dim]
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]

# # Save for Marcopolo
# sc <- import("scanpy", convert = FALSE)
# scvi <- import("scvi", convert = FALSE)
# library(anndata)
# raw_count = as.matrix(tiss[["RNA"]]@counts)[selected_genes,,drop = FALSE]
# tiss_adata <- sc$AnnData(
#   X = t(raw_count)
# )
# anndata::write_h5ad(tiss_adata, file.path(dir.path,data_source,"anndata_for_marcopolo",paste0(tissue_name,".h5ad")))

# Run Each Method and save results ========
dir.path <- dir.path0
folder.path = file.path(dir.path,"benchmark",data_source)
if(!file.exists(folder.path)){
  dir.create(folder.path, recursive = T)
}

# Create RunTime Table
if(!file.exists(file.path(folder.path, paste0(tissue_name,"_runtime.csv")))){
  df_runtime = data.frame("nGenes" = nrow(dat),
                          "nCells" = ncol(dat),
                          "Method" = method_ls,
                          'ncluster' = length(unique(tiss$annotation_group)),
                          row.names = method_ls)
}else{
  df_runtime = read.table(file.path(folder.path, paste0(tissue_name,"_runtime.csv")))
}
file.create(file.path(folder.path, "runtime_tmp.txt"))
for(method in method_ls){
  dir.file = file.path(folder.path, paste0(method,"_",tissue_name,".csv"))
  if(file.exists(dir.file)){
    next;
  }
  cat("Start Run ",method,"\n")
  start.time <- Sys.time()
  if(method == "lmd"){
    RunLMD(dat, feature_space, dir.file)
  }
  if(method == "hvg"){
    RunHVG(tiss, selected_genes, dir.file)
  }
  if(method == "wilcox_no_filter"){
    RunSeuratv4(tiss, selected_genes, feature_space, dir.file, res = 5)
  }
  if(method == "cosg"){
    RunCOSG(tiss, selected_genes, feature_space, dir.file, res = 5)
  }
  if(method == "scmarker"){
    RunSCMarker(dat, feature_space, dir.file)
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
  if(method == "marcopolo"){
    RunMarcopolo(tiss, selected_genes, dir.file)
  }
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  if(method %in% rownames(df_runtime)){
    df_runtime[method,"Time"] = as.numeric(time.taken, units = "mins")
  }else{
    tmp = df_runtime[1,]
    rownames(tmp) = tmp$Method = method
    tmp$'Time' = as.numeric(time.taken, units = "mins")
    df_runtime = rbind(df_runtime,tmp)
  }
}
file.remove(file.path(folder.path, "runtime_tmp.txt"))
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
if(!file.exists(folder.path)){
  dir.create(folder.path, recursive=T)
}
file_name = paste0(tissue_name,"_ground_truth_c1.txt")
if(!file.exists(file.path(folder.path,file_name))){
  avg_exp = AverageExpression(subset(tiss,features = selected_genes), assays = "RNA", slot = "counts", group.by = "annotation_group") %>% as.data.frame()
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
})

# # Concordance table
# celldb_marker = read.table(file.path(folder.path.gt, paste0(tissue_name,"_ground_truth_c2.txt")))[,1]
# df_benchmark = read.table(file.path(folder.path.rank, paste0(tissue_name,"_benchmark_rank_table.csv")),
#                           header = TRUE)
# # Among top N genes, how many of them are ground-truth?
# method_ls = c("lmd","hvg", "wilcox_no_filter",
#               "haystack","hotspot","semi","marcopolo")
# unlist(lapply(method_ls,function(method){
#   vec = setNames(df_benchmark[,method],df_benchmark$gene)
#   n = round(nrow(df_benchmark) * 0.05)
#   sum(names(sort(vec))[1:n] %in% celldb_marker)
# }))

# Cluster-based methods' adaptive resolution -------
method_subls = c("wilcox_no_filter","cosg")
lapply(tissue_ls, function(tissue_name){
  tiss = data_S_ls[[tissue_name]]
  n_dim = 20
  feature_space = as.matrix(tiss@reductions$pca@cell.embeddings[,1:n_dim])
  dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
  Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
  selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
  selected_genes = names(selected_genes)[selected_genes]
  dat = dat[selected_genes,,drop = FALSE]
  
  lapply(method_subls,function(method){
    folder.path = file.path(dir.path0,"benchmark",data_source,method)
    if(!file.exists(folder.path)){dir.create(folder.path)}
    if(!file.exists(file.path(folder.path, paste0(tissue_name,"_auroc.csv")))){
      res_grid = c(seq(0.5,3,0.5),seq(4,10,1))
      if(!file.exists(file.path(dir.path0,"benchmark",data_source,paste0(tissue_name,"_louvain_clusters.csv")))){
        df_res = unlist(lapply(res_grid,function(res){
          tiss_tmp = tiss
          tiss_tmp[["pca"]]@cell.embeddings = feature_space
          tiss_tmp <- tiss_tmp %>% FindNeighbors(dims = 1:n_dim, reduction = "pca") %>% 
            FindClusters(resolution = res, algorithm = 3)
          nlevels(Idents(tiss_tmp))
        }))
        df_res = data.frame(res = res_grid,
                            ncluster = df_res,
                            ncelltype = length(unique(tiss$'annotation_group')))
        write.table(df_res,file = file.path(dir.path0,"benchmark",data_source,paste0(tissue_name,"_louvain_clusters.csv")),row.names = FALSE)
        
      }
      df_res = read.table(file = file.path(dir.path0,"benchmark",data_source,paste0(tissue_name,"_louvain_clusters.csv")),header = T)
      res_grid = df_res$res[df_res$ncluster - df_res$ncelltype < 50]
      df_benchmark = lapply(res_grid, function(res){
        if(method == "wilcox_no_filter"){
          marker = RunSeuratv4(tiss, selected_genes, feature_space, res = res)
        }else if(method == "cosg"){
          marker = RunCOSG(tiss, selected_genes, feature_space, res = res)
        }
        if(!is.null(marker)){
          marker %>% select(gene,rank) %>% distinct()
        }else{
          marker
        }
      })
      names(df_benchmark) = res_grid
      df_benchmark = bind_rows(df_benchmark, .id = "res")
      df_benchmark <- tidyr::pivot_wider(df_benchmark, names_from = res, values_from = rank) %>% as.data.frame()
      df_benchmark[is.na(df_benchmark)] = nrow(df_benchmark)
      write.table(df_benchmark,file = file.path(folder.path, paste0(tissue_name,"_benchmark_rank_table.csv")),row.names = FALSE)
      # AUROC
      # df_benchmark = read.table(file.path(folder.path, paste0(tissue_name,"_benchmark_rank_table.csv")),
      #                           header = TRUE)
      folder.path.gt <- file.path(dir.path,data_source,"ground_truth_geneset")
      auc_df = do.call(rbind,lapply(colnames(df_benchmark)[-1],function(res){
        vec = quick_marker_benchmark(setNames(df_benchmark[,res],df_benchmark$gene),
                                     folder_path = folder.path.gt, 
                                     tissue_name = tissue_name)
        data.frame(gt_set = names(vec),AUROC = vec, res = res)
      }))
      write.table(auc_df,file = file.path(folder.path, paste0(tissue_name,"_auroc.csv")),row.names = FALSE)
    }
  })
})


# Density Index ========
dir.path = dir.path0
lapply(tissue_ls[c(1,5)],function(tissue_name){
tiss <- data_S_ls[[tissue_name]]
folder.path.rank <- file.path(dir.path,"benchmark",data_source)
df_benchmark = read.table(file.path(folder.path.rank, paste0(tissue_name,"_benchmark_rank_table.csv")),
                          header = TRUE)
rownames(df_benchmark) = df_benchmark$gene
n_dim = 20
low_qc = 1 # remove low-quality cells: cells express in less than N genes

# topn = c(seq(50,90,10),seq(100,900,100),seq(1000,3000,500),seq(4000,9000,1000),seq(10000,nrow(df_benchmark),2000))
topn = c(c(20,50),seq(100,1000,100),seq(1000,3000,500),seq(4000,nrow(df_benchmark),1000))
df_density_index =do.call(rbind,lapply(c(method_ls,"All genes"),function(method){
  if(method == "All genes"){
    top_gs = rownames(df_benchmark)
    topn = length(top_gs)
  }else if(method == "scmarker"){
    top_gs = rownames(df_benchmark)[df_benchmark[,method] == 1]
    topn = length(top_gs)
  }else{
    top_gs = rownames(df_benchmark)[order(df_benchmark[,method])]
  }
  result = do.call(rbind, lapply(topn, function(top){
    features = top_gs[1:top]
    # remove low-quality cells: cells express in less than N genes
    tmp_dat = subset(tiss,features = features)[["RNA"]]@counts
    pass_QC_cell = names(which(colSums(tmp_dat > 0) >= low_qc))
    subtiss = subset(tiss, cells = pass_QC_cell)
    subtiss <- subtiss %>% NormalizeData() %>% 
      ScaleData(features = features) %>% RunPCA(features = features, npcs = n_dim)
    # subtiss <- RunUMAP(subtiss, dims = 1:n_dim)
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
# df_density_index_null = do.call(rbind,lapply(topn,function(top){
#   set.seed(top)
#   seed_ls = sample(1:1e5,size = 100)
#   random_top_gs = lapply(seed_ls,function(sed){
#     set.seed(sed)
#     sample(rownames(df_benchmark),top,replace = FALSE)
#   })
#   result = unlist(lapply(random_top_gs, function(g_ls){
#     features = g_ls
#     # remove low-quality cells: cells express in less than N genes
#     tmp_dat = subset(tiss,features = features)[["RNA"]]@counts
#     pass_QC_cell = names(which(colSums(tmp_dat > 0) >= low_qc))
#     subtiss = subset(tiss, cells = pass_QC_cell)
#     
#     subtiss <- subtiss %>% ScaleData(features = features) %>% RunPCA(features = features, npcs = n_dim)
#     # Density-index
#     nn.dist = FNN::get.knn(Embeddings(subtiss, reduction = "pca"), k = 10, algorithm = "kd_tree")$nn.dist
#     sdVec <- stats::na.omit(subtiss@reductions$pca@stdev)
#     length_scale <- sqrt(sum(sdVec^2))
#     s2 = length_scale/mean(nn.dist)
#     s2
#   }))
#   result = data.frame(TopGene = top,DensityIndex = result)
# }) )
# write.table(df_density_index_null,file = file.path(folder.path.rank, paste0(tissue_name,"_DI_null.csv")),row.names = FALSE)

})



# Loop over LMD gene modules -------
tissue_name = "marrow_facs"
gene_module = lapply(tissue_name, function(tissue_name){
  local_gene = readRDS(file.path(folder.path,"intermediate_data",tissue_name,sprintf("local_gene_%s.rds",tissue_name)))
  local_gene$'gene_partition'
})[[1]]

df_benchmark = read.table(file.path(dir.path0,"benchmark",data_source, paste0(tissue_name,"_benchmark_rank_table.csv")),
                          header = TRUE)

pl = lapply(levels(gene_module),function(i){
  module = gene_module
  genes = names(module)[module == i]
  df = df_benchmark %>% filter(gene %in% genes)
  df = reshape2::melt(df,id = "gene",variable.name = "method",value.name = "rank")
  
  p = ggplot(df, aes(x = method, y = rank, fill = method)) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste0("Module ",i),
         x = "Method",
         y = "Rank") +
    theme(legend.position = "none")
  
})
p = wrap_plots(pl,nrow = 7)
ggsave("1.png",p,height = 28, width = 15)

