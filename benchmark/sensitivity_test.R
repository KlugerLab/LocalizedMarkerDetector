cran_packages <- c("readxl","dplyr","BiocManager")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)

# Load Data =======
dir.path <- "/banach1/ruiqi/local_marker/"
consortium = "tabular_muris"
tissue_name = "marrow_facs"
tiss <- readRDS(file.path(dir.path,"LMD_data",consortium,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]


# kNN sensitivity Test --------
folder.path = file.path(dir.path,"LMD_data",consortium,"intermediate_data",
                        tissue_name, "sensitivity_test")
dir.create(folder.path,recursive = T)
source(file.path(dir.path,"LocalMarkerDetector","LMD_function.R"))
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
write.table(score_df,file = file.path(folder.path,paste0(tissue_name,"_different_knn.csv")),row.names = TRUE)

# PC sensitivity Test --------
folder.path = file.path(dir.path,"LMD_data",consortium,"intermediate_data",
                        tissue_name, "sensitivity_test")
source(file.path(dir.path,"LocalMarkerDetector","LMD_function.R"))
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
write.table(score_df,file = file.path(folder.path,paste0(tissue_name,"_different_pc.csv")),row.names = TRUE)

# AUROC -----------
test_para_ls = c("pc","knn")
test_para = test_para_ls[2]
score_df = read.table(file.path(folder.path,paste0(tissue_name,sprintf("_different_%s.csv",test_para))))
score_ranked <- data.frame(score_df) %>% mutate(across(everything(), rank))
folder.path.gt <- file.path(dir.path,"LMD_data",consortium,"ground_truth_geneset")
auc_df = do.call(rbind,lapply(colnames(score_ranked),function(method){
  vec = quick_marker_benchmark(setNames(score_ranked[,method],rownames(score_ranked)),
                               folder_path = folder.path.gt, 
                               tissue_name = tissue_name)
  data.frame(gt_set = names(vec),AUROC = vec, Parameter = method)
}))
write.table(auc_df,file = file.path(folder.path, sprintf("%s_auroc_%s.csv",tissue_name,test_para)),row.names = FALSE)

# Jaccard Index
base = "ori_20PC"
base = "kNN5"
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
write.table(df,file = file.path(folder.path, sprintf("%s_jaccard_ind_%s.csv",tissue_name,test_para)),row.names = FALSE)

