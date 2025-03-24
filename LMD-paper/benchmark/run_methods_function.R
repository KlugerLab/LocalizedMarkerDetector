# Install Marcopolo: pip install marcopolo-pytorch --upgrade
# Install SEMITONES: pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip
# Install Hotspot: pip install git+https://github.com/yoseflab/Hotspot.git
# Install singleCellHaystack: install.packages("singleCellHaystack")

cran_packages <- c("reticulate","singleCellHaystack")
github_packages <- c("immunogenomics/presto","KChen-lab/SCMarker")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(github_packages, function(pkg) if (!requireNamespace(strsplit(pkg, "/")[[1]][2], quietly = TRUE)) remotes::install_github(pkg))
if (!requireNamespace("COSG", quietly = TRUE)) remotes::install_github("genecell/COSGR")
lapply(c(cran_packages, "presto","COSG","SCMarker"), require, character.only = TRUE)


#' LMD: Genes were ranked by the increasing order of Cumulative Diffusion-KL score
library(LocalizedMarkerDetector)
RunLMD <- function(dat, feature_space, dir.file){
  res = LMD(dat,feature_space,max_time = 2^20,knn = 5)
  marker = show_result_lmd(res)$'gene_table'
  write.table(marker,file = dir.file,row.names = TRUE)
}

#' HVG
RunHVG <- function(srat_obj, input_genes, dir.file){
  DefaultAssay(srat_obj) <- "RNA"
  marker <- subset(srat_obj,features = input_genes) %>% ScaleData(features = input_genes) %>% 
    FindVariableFeatures(nfeatures = length(input_genes)) %>% VariableFeatures() %>% as.data.frame() %>% setNames(c("gene")) %>% mutate(rank = 1:n())
  write.table(marker,file = dir.file,row.names = FALSE)
}

#' SinglecellHaystack: 
#' https://alexisvdb.github.io/singleCellHaystack/articles/examples/a02_example_scRNAseq.html
#' Genes were ranked based on the increasing order of log(adjusted P-value)
RunHaystack <- function(dat, feature_space, dir.file = NULL){
  set.seed(123)
  res <- haystack(x = feature_space, expression = dat)
  marker <- show_result_haystack(res.haystack = res)
  marker[order(marker$log.p.adj),"rank"] = 1:nrow(marker)
  if(is.null(dir.file)) return(marker)
  write.table(marker,file = dir.file,row.names = TRUE)
}

#' SEMITONES
#' https://github.com/ohlerlab/SEMITONES/tree/master
#' https://github.com/ohlerlab/SEMITONES_paper/blob/main/notebooks/10_SEMITONES_benchmark_corrected.ipynb
#' Genes were ranked by their lowest adjusted P-value in any reference cell (increasing order), if two genes had the same p.adj, prioritize genes with larger absolute enrichment score in any reference cell.
RunSEMITONES <- function(dat, feature_space, dir.file){
  cacheEnv <- new.env()
  options("reticulate.engine.environment" = cacheEnv)
  assign("dat",dat,envir = cacheEnv)
  assign("feature_space",feature_space,envir = cacheEnv)
  
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
  write.table(marker,file = dir.file, row.names = TRUE)
}

#' ## Hotspot
#' https://hotspot.readthedocs.io/en/latest/CD4_Tutorial.html
#' Genes were ranked by the increasing order FDR value (if two genes had the same, prioritize genes with larger Z-score).
RunHotspot <- function(srat_obj, input_genes, feature_space, dir.file = NULL){
  raw_count = as.matrix(srat_obj[["RNA"]]@counts)[input_genes,,drop = FALSE]
  cell_id = colnames(raw_count)
  gene_id = rownames(raw_count)
  total_count = srat_obj@meta.data$nCount_RNA
  
  cacheEnv <- new.env()
  options("reticulate.engine.environment" = cacheEnv)
  assign("raw_count",raw_count,envir = cacheEnv)
  assign("cell_id",cell_id,envir = cacheEnv)
  assign("gene_id",gene_id,envir = cacheEnv)
  assign("total_count",total_count,envir = cacheEnv)
  
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
  if(is.null(dir.file)) return(marker)
  write.table(marker,file = dir.file,row.names = TRUE)
}

#' ## SCMarker
#' https://github.com/KChen-lab/SCMarker
#' No rank is provided, only provide a set of genes
RunSCMarker <- function(dat, feature_space, dir.file = NULL){
  res = ModalFilter(data=dat,geneK=10,cellK=10,width=2)
  res = GeneFilter(res)
  res = getMarker(res, k = 300, n = 30)
  marker = data.frame(gene = rownames(dat), row.names = rownames(dat))
  marker[res$marker,"rank"] = 1
  marker$rank[is.na(marker$rank)] = nrow(marker)
  
  if(is.null(dir.file)){return(marker)}else{
    write.table(marker,file = dir.file,row.names = TRUE)
  }
}

#' ## COSG
#' https://github.com/genecell/COSGR/blob/main/vignettes/quick_start.Rmd
#' Genes were ranked by their largest COSG score in any cluster (decreasing order)
RunCOSG <- function(srat_obj, input_genes, feature_space, dir.file = NULL, res = 1.5){
  n_dim = dim(feature_space)[2]
  srat_obj[["pca"]]@cell.embeddings = feature_space
  srat_obj <- FindNeighbors(srat_obj, dims = 1:n_dim, reduction = "pca")
  srat_obj <- FindClusters(srat_obj, resolution = res, algorithm = 3)
  # if(nlevels(Idents(srat_obj)) > 50){
  #   return(NULL)
  # }
  marker <- cosg(
    subset(srat_obj,features = input_genes),
    groups='all',
    assay='RNA',
    slot='data',
    mu=1,
    n_genes_user=length(input_genes))
  marker <- dplyr::bind_cols(
    tidyr::pivot_longer(
      marker$names,
      cols = everything(),
      names_to = "cluster",
      values_to = "gene"
    ),
    dplyr::select(
      tidyr::pivot_longer(
        marker$scores,
        cols = everything(),
        names_to = "cluster",
        values_to = "raw_statistic"
      ),
      -cluster)
  )
  marker = marker  %>%
    group_by(gene) %>%
    filter(raw_statistic == max(raw_statistic)) %>%
    select(gene,raw_statistic) %>% distinct()
  marker[order(marker$raw_statistic,decreasing = TRUE),"rank"] = 1:nrow(marker)

  if(!is.null(dir.file)){
    write.table(marker,file = dir.file, row.names = TRUE)
  }else{
    return(marker)
  }
}

#' Seurat v4 Wilcoxon unfiltered: Genes were ranked by their lowest adjusted P-value in any cluster (increasing order), if two genes had the same, prioritize genes with larger FC.
RunSeuratv4 <- function(srat_obj, input_genes, feature_space, dir.file = NULL, res = 1.5){
  n_dim = dim(feature_space)[2]
  srat_obj[["pca"]]@cell.embeddings = feature_space
  srat_obj <- FindNeighbors(srat_obj, dims = 1:n_dim, reduction = "pca")
  srat_obj <- FindClusters(srat_obj, resolution = res, algorithm = 3)
  # if(nlevels(Idents(srat_obj)) > 50){
  #   return(NULL)
  # }
  marker <- presto:::wilcoxauc.Seurat(X = subset(srat_obj,features = input_genes),
                                      assay = "data", 
                                      seurat_assay = "RNA")
  colnames(marker)[1] = "gene"
  gene_rank = marker %>% arrange(.,padj,desc(logFC)) %>% distinct(.,gene, .keep_all = TRUE) %>% select(gene)
  gene_rank$'rank' = 1:nrow(gene_rank)
  marker <- merge(marker,gene_rank,by = "gene")
  if(!is.null(dir.file)){
    write.table(marker,file = dir.file, row.names = TRUE)
  }else{
    return(marker)
  }
}

#' Marcopolo
#' The rank of genes is given by `Marcopolo_rank`
RunMarcopolo <- function(srat_obj, input_genes, dir.file = NULL){
  raw_count = as.matrix(srat_obj[["RNA"]]@counts)[input_genes,,drop = FALSE]
  cacheEnv <- new.env()
  options("reticulate.engine.environment" = cacheEnv)
  assign("raw_count",raw_count,envir = cacheEnv)
  
  reticulate::py_run_string("
import pickle
import numpy as np
import pandas as pd
import torch
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import MarcoPolo

adata = sc.AnnData(X = r.raw_count.T)

if 'size_factor' not in adata.obs.columns:
    norm_factor = sc.pp.normalize_total(adata, exclude_highly_expressed=True, max_fraction= 0.2, inplace=False)['norm_factor']
    adata.obs['size_factor'] = norm_factor/norm_factor.mean()
    print('size factor was calculated')
    
regression_result = MarcoPolo.run_regression(adata=adata, size_factor_key='size_factor',
                         num_threads=1, device='cpu')

result = MarcoPolo.find_markers(adata=adata, regression_result=regression_result)
result = result.set_axis(list(adata.var.index) ,axis = 0)                           
                            ")
  marker = reticulate::py$result
  marker$'rank' = marker$MarcoPolo_rank + 1
  if(is.null(dir.file)) return(marker)
  write.table(marker,file = dir.file,row.names = TRUE)
}

RunSeuratv4_filter <- function(srat_obj, input_genes, feature_space, dir.file){
  srat_obj <- FindNeighbors(srat_obj, dims = 1:n_dim)
  srat_obj <- FindClusters(srat_obj, resolution = 1.5, algorithm = 3)
  marker <- FindAllMarkers(srat_obj, features = input_genes, only.pos = FALSE, test.use = "wilcox", logfc.threshold = 0.25,min.pct = 0.1)
  gene_rank = marker %>% arrange(.,p_val_adj,desc(avg_log2FC)) %>% distinct(.,gene, .keep_all = TRUE) %>% select(gene)
  gene_rank$'rank' = 1:nrow(gene_rank)
  marker <- merge(marker,gene_rank,by = "gene")
  write.table(marker,file = dir.file, row.names = TRUE)
}

quick_marker_benchmark <- function(gene_rank_vec,
                                   folder_path,
                                   tissue_name){
  max_logfc = read.table(file.path(folder_path, paste0(tissue_name,"_ground_truth_c1.txt"))) %>% rownames()
  celldb_marker = read.table(file.path(folder_path, paste0(tissue_name,"_ground_truth_c2.txt")))[,1]
  gt_list = c(lapply(seq(50,1000,50),function(x) max_logfc[1:x]),list(celldb_marker))
  names(gt_list) = c(paste0("Top",seq(50,1000,50)),"CellMarkerDB")
  df_benchmark = data.frame(gene_rank_vec,row.names = names(gene_rank_vec))
  
  auc_vec = do.call(c,lapply(1:length(gt_list),function(i){
    true_marker = gt_list[[i]]
    df_benchmark$"gt" = 0
    df_benchmark[true_marker,"gt"] = 1
    
    library(pROC)
    roc = roc(df_benchmark$gt, df_benchmark[,1], direction = ">")
    as.numeric(auc(roc))
  }))
  names(auc_vec) = names(gt_list)
  return(auc_vec)
}

quick_marker_benchmark_plot <- function(gene_rank_vec,
                                   folder_path,
                                   tissue_name){
  max_logfc = read.table(file.path(folder_path, paste0(tissue_name,"_ground_truth_c1.txt"))) %>% rownames()
  celldb_marker = read.table(file.path(folder_path, paste0(tissue_name,"_ground_truth_c2.txt")))[,1]
  gt_list = c(lapply(seq(50,1000,50),function(x) max_logfc[1:x]),list(celldb_marker))
  names(gt_list) = c(paste0("Top",seq(50,1000,50)),"CellMarkerDB")
  df_benchmark = data.frame(gene_rank_vec,row.names = names(gene_rank_vec))
  
  pl_roc = lapply(names(gt_list),function(i){
    true_marker = gt_list[[i]]
    df_benchmark$"gt" = 0
    df_benchmark[true_marker,"gt"] = 1
    
    library(pROC)
    roc = roc(df_benchmark$gt, df_benchmark[,1], direction = ">")
    
    roc_df <- data.frame(
      FPR = 1 - roc$specificities,  # False Positive Rate
      TPR = roc$sensitivities       # True Positive Rate
    )
    
    # Plot AUROC Curve
    p = ggplot(roc_df, aes(x = FPR, y = TPR)) +
      geom_line(color = "blue", size = 1) + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal reference line
      labs(title = paste0(i," (AUC = ", round(auc(roc), 3), ")"),
           x = "False Positive Rate (1 - Specificity)",
           y = "True Positive Rate (Sensitivity)") +
      theme_minimal()
    p
  })
  names(pl_roc) = names(gt_list)
  return(list(figure = pl_roc,
              gene_list = gt_list))
}
