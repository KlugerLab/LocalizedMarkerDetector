# Install Marcopolo: pip install marcopolo-pytorch --upgrade
# Install SEMITONES: pip install https://github.com/ohlerlab/SEMITONES/archive/master.zip
# Install Hotspot: pip install git+https://github.com/yoseflab/Hotspot.git
# Install singleCellHaystack: install.packages("singleCellHaystack")

cran_packages <- c("reticulate","singleCellHaystack")
github_packages <- c("immunogenomics/presto")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
# sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
sapply(github_packages, function(pkg) if (!requireNamespace(strsplit(pkg, "/")[[1]][2], quietly = TRUE)) remotes::install_github(pkg))
lapply(c(cran_packages), require, character.only = TRUE)

#' LMD: Genes were ranked by the increasing order of Cumulative Diffusion-KL score
dir.path <- "/banach1/ruiqi/local_marker/LocalMarkerDetector"
source(file.path(dir.path,"LMD_function.R"))
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

#' Seurat v4 Wilcoxon unfiltered: Genes were ranked by their lowest adjusted P-value in any cluster (increasing order), if two genes had the same, prioritize genes with larger FC.
RunSeuratv4 <- function(srat_obj, input_genes, feature_space, dir.file){
  n_dim = dim(feature_space)[2]
  srat_obj <- FindNeighbors(srat_obj, dims = 1:n_dim)
  srat_obj <- FindClusters(srat_obj, resolution = 1.5, algorithm = 3)
  marker <- presto:::wilcoxauc.Seurat(X = subset(srat_obj,features = input_genes),
                                      assay = "data", 
                                      seurat_assay = "RNA")
  colnames(marker)[1] = "gene"
  gene_rank = marker %>% arrange(.,padj,desc(logFC)) %>% distinct(.,gene, .keep_all = TRUE) %>% select(gene)
  gene_rank$'rank' = 1:nrow(gene_rank)
  marker <- merge(marker,gene_rank,by = "gene")
  write.table(marker,file = dir.file, row.names = TRUE)
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
