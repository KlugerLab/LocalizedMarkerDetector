cran_packages <- c("readxl","dplyr","BiocManager")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
if (!requireNamespace("GeneTrajectory", quietly = TRUE)) devtools::install_github("KlugerLab/GeneTrajectory")
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)
library(LocalizedMarkerDetector)
library(GeneTrajectory)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"run_methods_function.R"))
setwd("/banach1/ruiqi/local_marker/LocalizedMarkerDetector")
devtools::load_all(".") 

# mask the LocalizedMarkerDetector/GeneTrajectory function for testing purpose
coarse.grain <- function (cell.embedding, gene.expression, graph.dist, N = 1000, 
                          random.seed = 1) 
{
  message("Run k-means clustering")
  set.seed(random.seed)
  km.res <- stats::kmeans(cell.embedding, N, iter.max = 50)
  message("Coarse-grain matrices")
  KNN.membership.mat <- matrix(0, nrow = N, ncol = nrow(cell.embedding))
  for (i in 1:ncol(KNN.membership.mat)) {
    KNN.membership.mat[km.res$cluster[i], i] <- 1
  }
  KNN.membership.mat <- KNN.membership.mat/apply(KNN.membership.mat, 
                                                 1, sum)
  gene.expression.updated <- gene.expression %*% t(biclust::binarize(KNN.membership.mat, 0))
  graph.dist.updated <- KNN.membership.mat %*% graph.dist %*%
    t(KNN.membership.mat)
  
  colnames(graph.dist.updated) = rownames(graph.dist.updated) = paste0("cg",1:ncol(graph.dist.updated))
  colnames(gene.expression.updated) = colnames(graph.dist.updated)
  output <- list()
  output[["gene.expression"]] <- gene.expression.updated
  output[["graph"]] <- graph.dist.updated
  output
  return(output)
}

# Load Data =======
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")

dir.path <- dir.path0
data_source = "azimuth"
# tissue_name = "human_bone_marrow"
tissue_name = "human_lung_hlca"

folder.path = file.path(dir.path,data_source,"intermediate_data",
                        tissue_name, "sensitivity_test")
if(!dir.exists(folder.path)){
  dir.create(folder.path,recursive = T)
}

tiss <- readRDS(file.path(dir.path,data_source,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]
feature_space = Embeddings(tiss[["pca"]])[,1:20]
cell_graph = ConstructKnnGraph(feature_space = feature_space)

# Run Coarse graining -------
N.list = c(500, 1000, 3000, 5000, 10000, ncol(dat))
if(!file.exists(file.path(folder.path, paste0(tissue_name,"_cg_runtime.csv")))){
  df_runtime = data.frame("nCells" = N.list,
                          "Time_min" = NA)
}else{
  df_runtime = read.table(file.path(folder.path, paste0(tissue_name,"_cg_runtime.csv")))
}
file.create(file.path(folder.path, "runtime_tmp.txt"))

score_df = do.call(cbind,lapply(N.list, function(N){
  start.time <- Sys.time()
  
  if(N == ncol(dat)){
    W = cell_graph$graph
    dat_tmp = dat
  }else{
    cg_output = coarse.grain(feature_space,dat,cell_graph$graph,N = N)
    W = cg_output$graph
    dat_tmp = cg_output$gene.expression
  }
  
  rho = RowwiseNormalize(dat_tmp)
  rho = rho[, colnames(W), drop = FALSE]
  rho = rho[which(apply(rho, 1, function(x) sum(x > 0) >= 5)), 
            , drop = FALSE]
  res = fast_get_lmds(W = as.matrix(W), max_time = 2^20, init_state = rho)
  score_vec = res$cumulative_score
  score_vec[setdiff(selected_genes,names(score_vec))] = NA
  df = data.frame(score_vec[selected_genes])
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(paste0(N," ",as.numeric(time.taken, units = "mins")
             ,"mins,\n"), file = file.path(folder.path, "runtime_tmp.txt"), append = TRUE)
  
  # df_runtime[which(df_runtime$nCells == N),"Time_min"] = as.numeric(time.taken, units = "mins")
  
  df
}))

colnames(score_df) = N.list
write.table(score_df,file = file.path(folder.path,"different_cg.csv"),row.names = TRUE)
# write.table(df_runtime, file = file.path(folder.path,"runtime_cg.csv"),row.names = FALSE)
# file.remove(file.path(folder.path, "runtime_tmp.txt"))

score_df = read.table(file.path(folder.path,"different_cg.csv"), header = TRUE)
colnames(score_df) = gsub("X","",colnames(score_df))
score_ranked <- data.frame(score_df) %>% mutate(across(everything(), rank))
colnames(score_ranked) = gsub("X","",colnames(score_ranked))

## Jaccard Index ----------
base = as.character(ncol(dat))
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
colnames(df)[1:2] = c("Ncell","Top_Genes")
df$Top_Genes = as.factor(df$Top_Genes)
write.table(df,file = file.path(folder.path, "jaccard_ind_cg.csv"),row.names = FALSE)


## AUROC ----------
folder.path.gt <- file.path(dir.path,data_source,"ground_truth_geneset")
auc_df = do.call(rbind,lapply(colnames(score_ranked),function(method){
  vec = quick_marker_benchmark(setNames(score_ranked[,method],rownames(score_ranked)),
                               folder_path = folder.path.gt, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Parameter = method)
}))
write.table(auc_df,file = file.path(folder.path, "auroc_cg.csv"),row.names = FALSE)

