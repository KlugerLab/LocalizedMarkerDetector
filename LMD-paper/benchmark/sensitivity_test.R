cran_packages <- c("readxl","dplyr","BiocManager")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)
library(LocalizedMarkerDetector)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"run_methods_function.R"))

setwd("/banach1/ruiqi/local_marker/LocalizedMarkerDetector")
devtools::load_all(".")

# Load Data =======
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")

dir.path <- dir.path0
data_source = "tabular_muris"
tissue_name = "marrow_facs" # pancreas_facs, marrow_droplet
tiss <- readRDS(file.path(dir.path,data_source,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]

# Sensitivity Test ---------
folder.path = file.path(dir.path,data_source,"intermediate_data",
                        tissue_name, "sensitivity_test")
if(!dir.exists(folder.path)){
  dir.create(folder.path,recursive = T)
}
test_para_ls = c("pc","knn","cellprop","umi",
                 "binary","graph","ambient")
base_ls = c("ori_20PC","kNN5","cell_prop_1_1",
            "umi_prop_1_1","original","kNN","original")
names(base_ls) = test_para_ls

## kNN --------
test_para = test_para_ls[2]
# Use Default PCs
n_dim = 20
feature_space = Embeddings(tiss[["pca"]])[,1:n_dim]
score_df = do.call(cbind,lapply(c(3:10,20,30,40,50),function(kNN){
  try({
    res = LMD(dat,feature_space,knn = kNN)
    return(res$cumulative_score)
  },silent = TRUE)
}))
colnames(score_df) = paste0("kNN",c(3:10,20,30,40,50))
write.table(score_df,file = file.path(folder.path,sprintf("%s_different_%s.csv",tissue_name,test_para)),row.names = TRUE)


## PC  --------
test_para = test_para_ls[1]
# Generate PC embeddings
DefaultAssay(tiss) = "RNA"
tmp_tiss = tiss %>% RunPCA(npcs = 100, verbose = FALSE)
feature_space_all = Embeddings(tmp_tiss[["pca"]]); rm(tmp_tiss)
pc_list = c("ori_20PC","ori_tSNE",10,20,30,40,50,75,100)
score_df = do.call(cbind,lapply(pc_list,function(n_dim){
  try({
    if(n_dim == "ori_20PC"){
      feature_space = Embeddings(tiss[["pca"]])[,1:20]
    }else if(n_dim == "ori_tSNE"){
      feature_space = Embeddings(tiss[["tsne"]])[,1:2]
    }else{
      n_dim = as.numeric(n_dim)
      feature_space = feature_space_all[,1:n_dim]
    }
    res = LMD(dat, feature_space = feature_space, knn = 5)
    return(res$cumulative_score)
  },silent = TRUE)
}))
colnames(score_df) = c(pc_list[c(1,2)],paste0(pc_list[-c(1,2)],"PC"))
write.table(score_df,file = file.path(folder.path,sprintf("%s_different_%s.csv",tissue_name,test_para)),row.names = TRUE)

## Cell number with bootstrap -------
test_para = test_para_ls[3]
n_dim = 20 # default PCs
knn = 5 # default knn

cell_prop_ls = c(0.01,0.05,0.1,0.25,0.5,0.75,1)
nrep = 10 # set the # of subsampling for each cell number
seed_ls = lapply(1:length(cell_prop_ls),function(i){
  set.seed(i)
  sample(1:1e5,nrep,replace = FALSE)})
names(seed_ls) =  cell_prop_ls 
score_df = do.call(cbind,lapply(cell_prop_ls,function(cell_prop){
  try({
    ncell = floor(ncol(tiss) * cell_prop)
    nrep_to_run <- if (cell_prop == 1) 1 else 1:nrep
    df = do.call(cbind, lapply(nrep_to_run, function(i){
      set.seed(seed_ls[[as.character(cell_prop)]][i])
      cells = sample(colnames(tiss))[1:ncell]
      tmp_tiss = subset(tiss, cells = cells) %>% 
        RunPCA(npcs = n_dim, verbose = FALSE)
      feature_space_tmp = Embeddings(tmp_tiss[["pca"]])[,1:n_dim]
      dat_tmp = dat[,colnames(tmp_tiss)]
      res = LMD(dat_tmp,feature_space_tmp,knn = knn)
      score_vec = res$cumulative_score
      score_vec[setdiff(selected_genes,names(score_vec))] = NA
      df = data.frame(score_vec[selected_genes]); colnames(df) = paste0("rep",i)
      df
    }))
    colnames(df) = paste0("rep",nrep_to_run)
    df
  },silent = TRUE)
}))
colnames(score_df) = c(paste0(rep(paste0("cell_prop_",setdiff(cell_prop_ls,1)),each = nrep),"_",
                            rep(1:nrep,length(cell_prop_ls)-1)),"cell_prop_1_1")
write.table(score_df,file = file.path(folder.path,sprintf("%s_different_%s.csv",tissue_name,test_para)),row.names = TRUE)

## Data sparsity (UMIs) ------
#' https://github.com/const-ae/transformGamPoi-Paper/blob/master/benchmark/src/downsampling_benchmark/download_deeply_sequenced_datasets.R
test_para = test_para_ls[4]

n_dim = 20 # default PCs
knn = 5 # default knn
count_mtx = tiss[["RNA"]]@counts
  
umi_prop_ls = c(0.005, 0.01, 0.05, 0.1, 1, 0.0005, 0.001, 0.00015, 0.0001)
nrep = 10 # set the # of subsampling for each cell number
seed_ls = lapply(1:length(umi_prop_ls),function(i){
  set.seed(i)
  sample(1:1e5,nrep,replace = FALSE)})
names(seed_ls) =  umi_prop_ls 
score_df = do.call(cbind,lapply(umi_prop_ls[6:7],function(umi_prop){
  try({
    nrep_to_run <- if (umi_prop == 1) 1 else 1:nrep
    df = do.call(cbind, lapply(nrep_to_run, function(i){
      # Original data, downsample UMI
      set.seed(seed_ls[[as.character(umi_prop)]][i])
      downsampled_count_mtx = scuttle::downsampleMatrix(count_mtx, prop = umi_prop, bycol = TRUE)
      # Process the downsampled data
      tmp_tiss = CreateSeuratObject(downsampled_count_mtx)
      VariableFeatures(tmp_tiss) = VariableFeatures(tiss)
      tmp_tiss <- tmp_tiss %>% NormalizeData(scale.factor = 1e6) %>% ScaleData() %>%
        RunPCA(npcs = n_dim, verbose = FALSE)
      
      # Extract expression matrix & feature_space
      dat_tmp = as.matrix(tmp_tiss[["RNA"]]@data)
      dat_tmp = dat_tmp[rownames(dat),colnames(dat)]
      dat_tmp = dat_tmp[matrixStats::rowSums2(dat_tmp) > 0,
                        matrixStats::colSums2(dat_tmp) > 0]
      feature_space_tmp = Embeddings(tmp_tiss[["pca"]])[,1:n_dim]
      feature_space_tmp = feature_space_tmp[colnames(dat_tmp),]

      # RunLMD
      res = LMD(dat_tmp,feature_space_tmp,knn = knn)
      score_vec = res$cumulative_score
      score_vec[setdiff(selected_genes,names(score_vec))] = NA
      df = data.frame(score_vec[selected_genes]); colnames(df) = paste0("rep",i)
      df
    }))
    colnames(df) = paste0("rep",nrep_to_run)
    df
  },silent = TRUE)
}))
colnames(score_df) = c(paste0(rep(paste0("umi_prop_",setdiff(umi_prop_ls,1)),each = nrep),"_",
                              rep(1:nrep,length(umi_prop_ls)-1)),"umi_prop_1_1")
write.table(score_df,file = file.path(folder.path,sprintf("%s_different_%s.csv",tissue_name,test_para)),row.names = TRUE)

## Binary expression --------
test_para = test_para_ls[5]
n_dim = 20
feature_space = Embeddings(tiss[["pca"]])[,1:n_dim]
score_df = do.call(cbind,lapply(c("original","binary"), function(x){
  dat_tmp = dat
  if(x == "binary"){
    dat_tmp[dat_tmp > 0] = 1
  }
  res = LMD(dat_tmp,feature_space,knn = 5)
  score_vec = res$cumulative_score
  score_vec[setdiff(selected_genes,names(score_vec))] = NA
  df = data.frame(score_vec[selected_genes])
  df
}))
colnames(score_df) = c("original","binary")
write.table(score_df,file = file.path(folder.path,sprintf("%s_different_%s.csv",tissue_name,test_para)),row.names = TRUE)


## Alternative Graphs (load different datasets) --------
test_para = test_para_ls[6]

dir.path <- dir.path0
data_source = "tabular_muris"
tissue_ls = c("marrow_facs","pancreas_facs","lung_facs","marrow_droplet")
data_source = "azimuth"
tissue_ls = c("human_kidney","human_pancreas","mouse_motor_cortex")

lapply(tissue_ls,function(tissue_name){
  tiss <- readRDS(file.path(dir.path,data_source,paste0(tissue_name,".rds")))
  DefaultAssay(tiss) <- "RNA"
  dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
  Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
  selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
  selected_genes = names(selected_genes)[selected_genes]
  dat = dat[selected_genes,,drop = FALSE]
  folder.path = file.path(dir.path,data_source,"intermediate_data",
                          tissue_name, "sensitivity_test")
  if(!dir.exists(folder.path)){
    dir.create(folder.path,recursive = T)
  }
  n_dim = 20
  feature_space = Embeddings(tiss[["pca"]])[,1:n_dim]
  score_df = do.call(cbind,lapply(c("kNN","SNN","Gaussian"), function(method){
    if(method == "kNN"){
      res = LMD(dat,feature_space,knn = 5,kernel = FALSE)
    }else{
      res = LMD(dat,feature_space,knn = 5,kernel = TRUE, kernel_used = method)
    }
    score_vec = res$cumulative_score
    score_vec[setdiff(selected_genes,names(score_vec))] = NA
    df = data.frame(score_vec[selected_genes])
    df
  }))
  colnames(score_df) = c("kNN","SNN","Gaussian")
  write.table(score_df,file = file.path(folder.path,sprintf("%s_different_%s.csv",tissue_name,test_para)),row.names = TRUE)
  
})

## Ambient RNA --------
test_para = test_para_ls[7]
folder.path = file.path(dir.path,"simulated_pancreas", "sensitivity_test")
if(!dir.exists(folder.path)){
  dir.create(folder.path,recursive = T)
}
data_S_ls = list(original = readRDS(file.path(folder.path,"mislet_simulated.rds")),
                 low = readRDS(file.path(folder.path,"synthetic_contaminated_low.rds")),
                 mid = readRDS(file.path(folder.path,"synthetic_contaminated_medium.rds")),
                 high = readRDS(file.path(folder.path,"synthetic_contaminated_high.rds")))
disturbed_genes = readRDS(file.path(folder.path,"gene500.rds"))

# Define genes used for testing
data_S <- data_S_ls$ground_truth
DefaultAssay(data_S) <- "RNA"
dat = as.matrix(data_S[[DefaultAssay(data_S)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.8)
selected_genes = names(selected_genes)[selected_genes]
selected_genes = union(selected_genes,disturbed_genes)
rm(dat,data_S,Gene_detected_count)

score_df = do.call(cbind,lapply(data_S_ls, function(data_S){
  dat = as.matrix(data_S[[DefaultAssay(data_S)]]@data)[selected_genes,]
  n_dim = 10
  feature_space = Embeddings(data_S[["pca"]])[,1:n_dim]
  res = LMD(dat, feature_space, knn = 10, kernel = FALSE)
  score_vec = res$cumulative_score
  score_vec[setdiff(selected_genes,names(score_vec))] = NA
  df = data.frame(score_vec[selected_genes])
  df
}))
colnames(score_df) = names(data_S_ls)
write.table(score_df,file = file.path(folder.path,sprintf("different_%s.csv",test_para)),row.names = TRUE)
score_df = read.table(file.path(folder.path,sprintf("different_%s.csv",test_para)))
score_ranked <- data.frame(score_df) %>% mutate(across(everything(), rank))


# Metric -------
## AUROC -----------
# test_para = test_para_ls[5]
score_df = read.table(file.path(folder.path,paste0(tissue_name,sprintf("_different_%s.csv",test_para))))
score_ranked <- data.frame(score_df) %>% mutate(across(everything(), rank))

folder.path.gt <- file.path(dir.path,data_source,"ground_truth_geneset")
auc_df = do.call(rbind,lapply(colnames(score_ranked),function(method){
  vec = quick_marker_benchmark(setNames(score_ranked[,method],rownames(score_ranked)),
                               folder_path = folder.path.gt, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Parameter = method)
}))
write.table(auc_df,file = file.path(folder.path, sprintf("%s_auroc_%s.csv",tissue_name,test_para)),row.names = FALSE)

## Jaccard Index ----------
base = base_ls[test_para]
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
colnames(df)[1:2] = c(test_para,"Top_Genes")

if(test_para == "ambient"){
  write.table(df,file = file.path(folder.path, sprintf("jaccard_ind_%s.csv",test_para)),row.names = FALSE)
}else{
  write.table(df,file = file.path(folder.path, sprintf("%s_jaccard_ind_%s.csv",tissue_name,test_para)),row.names = FALSE)
}
