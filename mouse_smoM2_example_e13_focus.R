dir.path <- "/banach1/ruiqi/local_marker/LocalMarkerDetector"
source(file.path(dir.path,"LMD_function.R"))

# Define package lists
bioc_packages <- c("clusterProfiler", "AnnotationDbi", "ReactomePA", "org.Mm.eg.db",
                   "gprofiler2", "msigdbr")
github_packages <- c("satijalab/seurat-wrappers");
names(github_packages) = strsplit(github_packages, "/")[[1]][2]
names(github_packages)[grepl("seurat-wrappers",names(github_packages))] = "SeuratWrappers"
# Install packages
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
sapply(1:length(github_packages), function(i) if (!requireNamespace(names(github_packages)[i], quietly = TRUE)) remotes::install_github(github_packages[i]))
lapply(c(bioc_packages,names(github_packages)), require, character.only = TRUE)

# Load Data --------------------
sample_ls = unlist(lapply(c("E13.5", "E14.5"),function(i) paste0(i,"_", c("MUT", "CTL"))))
sample_ls = paste0("smom2_dermal_",sample_ls)
for(i in sample_ls){
  assign(paste0("data_S_",i), 
         readRDS(sprintf("/data/rq25/peggy_10x_mouse_smom2_202112/seurat_analysis/result/data/data_S_peggy_%s.rds",i)))
}
# Preprocess ------------
dir.path = "/banach1/ruiqi/Peggy_data/Peggy_scdata/240329_smom2"
folder.path <- file.path(dir.path,"process_data")
dir.create(folder.path, recursive=T)
lapply(sample_ls, function(sample_name){
  data_S <- get(paste0("data_S_",sample_name))
  data_S <- data_S %>% NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
  min.pc = FindPC(srat = data_S)
  data_S <- RunUMAP(data_S, dims = 1:min.pc, seed.use = 42)
  data_S <- data_S %>% FindNeighbors() %>% FindClusters(res = 1.5)
  assign(paste0("data_S_",sample_name),data_S, envir = .GlobalEnv)
})

# Trim & Re-embed data
#' Remove mosaics from E14.5 MUT 
sample_name = sample_ls[3]
data_S <- get(paste0("data_S_",sample_name))
data_S = subset(data_S, RNA_snn_res.1.5 != 16)
data_S <- data_S %>% NormalizeData() %>% FindVariableFeatures() %>%
  ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)
min.pc = FindPC(srat = data_S)
data_S <- RunUMAP(data_S, dims = 1:min.pc, seed.use = 42)
data_S <- data_S %>% FindNeighbors(dims = 1:min.pc) %>% 
  FindClusters(dims = 1:min.pc, res = 1.5)
assign(paste0("data_S_",sample_name),data_S)

# Add CellType Annotation -----------
#' Annotate based on Dkk1, Dkk2, Lef1, Sox2
sample_name = sample_ls[1]
data_S <- get(paste0("data_S_",sample_name))
p1 = DimPlot(data_S, group.by = "RNA_snn_res.1.5",label = TRUE)
p2 = FeaturePlot(data_S,label = TRUE, order = TRUE, features = c("Dkk1","Lef1","Dkk2","Sox2")) & NoAxes()
p1 + p2
# E13.5 MUT
data_S$'celltype' = ifelse(data_S$RNA_snn_res.1.5 %in% c(11,10,13,2,4,9), "LD", 
                           ifelse(data_S$RNA_snn_res.1.5 %in% 14,"DC","UD"))
# E13.5 CTL
data_S$'celltype' = ifelse(data_S$RNA_snn_res.1.5 %in% c(1,2,3,4,7,9,17,12,18,19), "LD", "UD")
# E14.5 MUT
data_S$'celltype' = ifelse(data_S$RNA_snn_res.1.5 %in% c(0,11,12), "LD", 
                           ifelse(data_S$RNA_snn_res.1.5 %in% c(13,15),"DC","UD"))
# E14.5 CTL
data_S$'celltype' = ifelse(data_S$RNA_snn_res.1.5 %in% c(0,2,8,11,13) & data_S@reductions$umap@cell.embeddings[,1] < 0, "LD", 
                           ifelse(data_S$RNA_snn_res.1.5 %in% 7,"DC","UD"))
DimPlot(data_S, group.by = "celltype")
data_S <- RunALRA(data_S); DefaultAssay(data_S) <- "RNA"
saveRDS(data_S,file = file.path(folder.path,sprintf("data_S_%s.rds",sample_name)))

# RunLMD ---------------
folder.path <- file.path(dir.path,"intermediate_data")
dir.create(folder.path, recursive=T)
lapply(sample_ls,function(sample_name){
  data_S <- get(sprintf("data_S_%s",sample_name))
  cat(sample_name,":",FindPC(data_S),"PC\n")
  dat = as.matrix(data_S[["RNA"]]@data)
  Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
  selected_genes = names(Gene_detected_count)[(Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)]
  dat = dat[selected_genes,,drop = FALSE]
  feature_space = Embeddings(data_S[["pca"]])[,1:FindPC(srat = data_S)]
  tic()
  lmd_result = LMD(dat, feature_space, max_time = 2^20)
  toc()
  local_gene = show_result_lmd(lmd_result)
  saveRDS(local_gene,file = file.path(folder.path,paste0("local_gene_",sample_name,".rds")))
})

# Group_Localized_genes ----------------
lapply(sample_ls,function(sample_name){
  data_S <- get(sprintf("data_S_%s",sample_name))
  local_gene <- readRDS(file.path(folder.path,paste0("local_gene_",sample_name,".rds")))
  dat_alra = as.matrix(data_S[["alra"]]@data)[names(local_gene$cut_off_gene),]
  dist = Calculate_distance(dat_alra, method = "jaccard")
  saveRDS(dist, file = file.path(folder.path,paste0(sample_name,"_dist.rds")))
  if(sample_name == sample_ls[1]){
    deepSplit = 2}else{deepSplit = 1}
  res = Obtain_gene_partition(dist, clustering_method = "average", 
                              deepSplit = deepSplit, return_tree = TRUE)
  saveRDS(res, file = file.path(folder.path,paste0(sample_name,"_gene_tree.rds")))
  # Visualize_gene_heatmap(dist, gene_partition = res$gene_partition, gene_hree = res$gene_hree)
  Visualize_gene_tsne(dist, gene_partition = res$gene_partition)
  pl = FeaturePlot_meta(data_S, feature_partition = res$gene_partition)
})
# manually merge some modules:
lapply(sample_ls[4],function(sample_name){
  res = readRDS(file = file.path(folder.path,paste0(sample_name,"_gene_tree.rds")))
  if(sample_name == sample_ls[4]){
    # merge wnt-modules
    levels(res$gene_partition)[levels(res$gene_partition)==6] = 5
    # merge dc-modules
    levels(res$gene_partition)[levels(res$gene_partition)==3] = 2
  }
  saveRDS(res, file.path(folder.path,paste0(sample_name,"_gene_tree.rds")))
})

# Calculate Module score ----------
gene_partition_all = unlist(lapply(sample_ls,function(sample_name){
  res = readRDS(file = file.path(folder.path,paste0(sample_name,"_gene_tree.rds")))
  gene_partition = res$gene_partition
  levels(gene_partition) = paste0(sample_name,"_",levels(gene_partition))
  return(gene_partition)
}))

cell_block_ls = lapply(sample_ls,function(sample_name){
  data_S <- get(sprintf("data_S_%s",sample_name))
  data_S@meta.data = data_S@meta.data[,!grepl("Module",colnames(data_S@meta.data))]
  data_S = AddModuleActivityScore(data_S, gene_partition = gene_partition_all, do_local_smooth = FALSE)
  # Obtain ModuleScore mtx
  col_id = grep("Module",colnames(data_S@meta.data))
  cell_block = data_S@meta.data[,col_id]
  cell_block
})
names(cell_block_ls) = sample_ls
saveRDS(cell_block_ls, file = file.path(folder.path,"cell_block_ls.rds"))

# Filter-out modules express in less than 10 cells
modules_rename_ls = unlist( lapply(sample_ls,function(sample_name){
  cell_block = cell_block_ls[[sample_name]]
  col_id = grep(sample_name,colnames(cell_block))
  cell_block = cell_block[,col_id]
  modules = names(which(colSums(cell_block > 0.5) >= 10))
  prefix = unique(sub("_(\\d+)$", "_", modules))
  modules = modules[order(as.numeric(sub(prefix, "", modules)))]
  modules_rename = setNames(paste0(prefix,1:length(modules)),modules)
  modules_rename
}) )

# Update saved partition results
cell_block_ls = lapply(cell_block_ls, function(cell_block){
  cell_block = cell_block[,colnames(cell_block) %in% names(modules_rename_ls)]
  colnames(cell_block) = modules_rename_ls[colnames(cell_block)]
  cell_block
})
lapply(sample_ls,function(sample_name){
  res = readRDS(file.path(folder.path,paste0(sample_name,"_gene_tree.rds")))
  local_gene = readRDS(file.path(folder.path,paste0("local_gene_",sample_name,".rds")))
  res$gene_partition_rename = res$gene_partition
  levels(res$gene_partition_rename) = 
    modules_rename_ls[paste0("Module",sample_name,"_",levels(res$gene_partition_rename))]
  levels(res$gene_partition_rename) = sub(".*_(\\d+)$", "\\1", levels(res$gene_partition_rename))
  res$gene_partition_rename = res$gene_partition_rename[!is.na(res$gene_partition_rename)]
  local_gene$'gene_partition' = res$gene_partition_rename
  saveRDS(res, file = file.path(folder.path,paste0(sample_name,"_gene_tree.rds")))
  saveRDS(local_gene, file = file.path(folder.path,paste0("local_gene_",sample_name,".rds")))
  
  data_S <- get(sprintf("data_S_%s",sample_name))
  data_S@meta.data = data_S@meta.data[,!grepl("Module",colnames(data_S@meta.data))]
  data_S <- AddMetaData(data_S, cell_block_ls[[sample_name]])
  assign(paste0("data_S_",sample_name), data_S, envir = .GlobalEnv)
  saveRDS(data_S,file = file.path(dir.path,"process_data",sprintf("data_S_%s.rds",sample_name)))
})

# Print genes of each module
lapply(sample_ls,function(sample_name){
  local_gene = readRDS(file.path(folder.path,paste0("local_gene_",sample_name,".rds")))
  cat(sample_name,"\n")
  for(i in levels(local_gene$gene_partition)){
    cat("Module",i,"\n")
    message(paste0(names(local_gene$gene_partition)[local_gene$gene_partition == i],collapse = ", "))
  }
})

# Calculate Jaccard Index of E13.5 mut modules with modules from other samples -------
modules = setNames(paste0("Module",sample_ls[1],"_",c(1,2,3,12,14,15,16)),
                   c(paste0("CC Module",1:3),"Quiescent Module","Wnt Module","DC Module", "High-Wnt Module"))
jaccard_index_ls = lapply(sample_ls[2:4],function(sample_name){
  data_S = get(paste0("data_S_",sample_name))
  col_names = colnames(data_S@meta.data)[grepl("Module",colnames(data_S@meta.data))]
  col_names = col_names[col_names %in% modules | grepl(sample_name,col_names)]
  jaccard_index <- 1 - philentropy::distance(
    t(data_S@meta.data[,col_names]), 
    method = "jaccard", mute.message=TRUE)
  rownames(jaccard_index) = colnames(jaccard_index) = col_names
  jaccard_index = jaccard_index[col_names %in% modules,
                                !col_names %in% modules,drop = FALSE]
  jaccard_index
}); names(jaccard_index_ls) = sample_ls[2:4]
df = do.call(rbind,lapply(names(jaccard_index_ls),function(sample_name){
  jaccard_index = jaccard_index_ls[[sample_name]]
  df = data.frame(module1 = rownames(jaccard_index),
             module2 = colnames(jaccard_index)[apply(jaccard_index,1,which.max)],
             jaccard_index = apply(jaccard_index,1,max),
             sample = sample_name)
  df %>% filter(jaccard_index > 0.35)
}))
write.csv(df,file = file.path(folder.path,"module_colocalize_align.csv"),row.names = FALSE)

# Rename some specific modules -------
df$'nickname' = names(modules)[match(df$module1,modules)]
module_ls = list(
  modules,
  setNames(paste0("Module",sample_ls[2],"_",c(17,14,16,15)),
           c("Wnt Module",paste0("CC Module",1:3))),
  setNames(paste0("Module",sample_ls[3],"_",c(15)),
           "Quiescent Module"),
  setNames(paste0("Module",sample_ls[4],"_",c(3,1)),
           c("Wnt Module","DC Module"))
); names(module_ls) = sample_ls[1:4]
saveRDS(module_ls, file = file.path(folder.path,"module_ls.rds"))

# Boxplot of ModuleActivityScore of active cells -------
modules = module_ls$smom2_dermal_E13.5_MUT[c("DC Module","High-Wnt Module","Wnt Module","Quiescent Module")]
nn_score_df = do.call(rbind,lapply(modules,function(module_name){
  nn_score_df = do.call(rbind,lapply(sample_ls,function(sample_name){
    data_S = get(paste0("data_S_",sample_name))
    if(!module_name %in% colnames(data_S@meta.data)){return(NULL)}
    df = do.call(rbind,lapply(c(5,10,15,20),function(knn){
      A = Symmetric_KNN_graph(knn = knn, Embeddings(data_S@reductions$pca)[,1:FindPC(data_S)])$adj_matrix
      diag(A) = 0
      a = ((A %*% data_S@meta.data[,module_name]) / rowSums(A))[data_S@meta.data[,module_name] > 0.5,]
      data.frame(value = a, sample = sample_name, knn = knn)
    }))
    df
  }))
  nn_score_df$sample = as.factor(nn_score_df$sample)
  nn_score_df$knn = as.factor(nn_score_df$knn)
  nn_score_df$'module' = module_name
  nn_score_df
}))
nn_score_df$'module' = names(modules)[match(nn_score_df$'module',modules)]
write.csv(nn_score_df, file = file.path(folder.path,"neighbors_module_score.csv"))

# Pathway Enrichment Analysis for each sample---------
library("ReactomePA")
library("org.Mm.eg.db")
library("clusterProfiler")
library("msigdbr")
pathway_result = lapply(sample_ls,function(sample_name){
  local_gene = get(paste0("local_gene_",sample_name))
  
  # Set all genes(after filtering) in the expression matrix as background
  selected_genes = names(local_gene$gene_rank)
  universe_df = data.frame("symbol" = selected_genes,"entrez" = mapIds(org.Mm.eg.db, keys=selected_genes, column="ENTREZID", keytype="SYMBOL"))
  universe_df = universe_df[!is.na(universe_df$entrez),]
  gene_partition = local_gene$gene_partition
  
  # Reactome Analysis
  epathway_result <- lapply(levels(gene_partition), function(i){
    enrichPathway(gene=universe_df[universe_df$symbol %in% names(gene_partition)[gene_partition == i],"entrez"],
                  organism = "mouse",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  universe = universe_df$entrez)
  })
  names(epathway_result) = levels(gene_partition)
  epathway_result <- Filter(Negate(is.null), epathway_result)
  
  # GO (MF) Analysis
  ego_result <- lapply(levels(gene_partition), function(i){
    enrichGO(gene = universe_df[universe_df$symbol %in% names(gene_partition)[gene_partition == i],"entrez"],
             OrgDb = 'org.Mm.eg.db', # mouse
             keyType = "ENTREZID",
             ont = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             universe = universe_df$entrez)
  })
  names(ego_result) = levels(gene_partition)
  ego_result <- Filter(Negate(is.null), ego_result)
  
  # Msig Analysis
  m_df = msigdbr(species = "Mus musculus", category = "H") %>% 
    dplyr::select(gs_name, entrez_gene)
  m_df$gs_name = gsub("HALLMARK_","",m_df$gs_name)
  msig_result <- lapply(levels(gene_partition), function(i){
    enricher(gene=universe_df[universe_df$symbol %in% names(gene_partition)[gene_partition == i],"entrez"],
             TERM2GENE = m_df,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             universe = universe_df$entrez)
  })
  names(msig_result) = levels(gene_partition)
  msig_result <- Filter(Negate(is.null), msig_result)
  
  ls = list(reactome = epathway_result,
       go = ego_result,
       msig = msig_result)
  return(ls)
})
names(pathway_result) = sample_ls
saveRDS(pathway_result, file = file.path(folder.path,"pathway_result.rds"))

# Regressing-out Cell Cycle embedding ------- 
lapply(sample_ls[c(1,3)],function(sample_name){
  data_S = get(paste0("data_S_",sample_name))
  data_S <- ScaleData(data_S, vars.to.regress = c("S.Score", "G2M.Score"), 
                      features = rownames(data_S))
  data_S <- data_S %>% RunPCA(features = VariableFeatures(data_S)) %>% RunUMAP(dims = 1:FindPC(srat = data_S), seed.use = 42)
  data_tmp = get(paste0("data_S_",sample_name))
  data_tmp@reductions[["umap_cc_regressout"]] = data_S@reductions$umap
  assign(paste0("data_S_",sample_name),data_tmp,envir = .GlobalEnv); rm(data_tmp)
})

# Merged Cell Cycle embedding -------  
library(gprofiler2)
s.genes = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
g2m.genes = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
genes_CC = c(s.genes,g2m.genes)

### Merge MUT & CTL
data_S = lapply(sample_ls[1:2],function(sample_name){
  get(paste0("data_S_",sample_name))
})
data_S = merge(data_S[[1]],data_S[-1])
data_S$'type' = ifelse(data_S$old.ident == "E13_SmoM2","E13.5_SmoM2","E13.5_WT")
data_S <- data_S %>% NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData(verbose = FALSE, features = rownames(data_S))
data_S <- data_S  %>% 
  RunPCA(npcs = 50, verbose = FALSE, features = VariableFeatures(data_S))
data_S <- data_S %>%
  RunUMAP(dims = 1:FindPC(srat = data_S), seed.use = 42, reduction.name = "umap_hvg")
data_S <- data_S  %>% 
  RunPCA(npcs = 50, verbose = FALSE, features = genes_CC)
data_S <- data_S %>%
  RunUMAP(dims = 1:FindPC(srat = data_S), seed.use = 42, reduction.name = "umap_cc")
data_S <- DietSeurat(
  data_S, assays = "RNA", dimreducs = c("umap_hvg","umap_cc")
)
saveRDS(data_S, file = file.path(dir.path,"process_data","data_S_smom2_dermal_E13.5.rds"))


# Wnt level -------
sample_name = sample_ls[1]
modules = module_ls[[sample_name]]
modules = modules[grepl("DC|Wnt",names(modules))]
modules = modules[c(1,3,2)]
data_S = get(paste0("data_S_",sample_name))
data_S$'wnt_related_module' = "Other"
thred = 0.4
for(m in names(modules)){
  data_S$'wnt_related_module' = ifelse(data_S@meta.data[,modules[m]] > thred, m, data_S$'wnt_related_module')
}
genes = "Lef1"
df = FetchData(data_S, vars = c(genes, "wnt_related_module"))
df$Phase = data_S$Phase; 
df = cbind(df,data_S@meta.data[,modules])
write.csv(df,file = file.path(folder.path,"e13_smom2_wnt_level_table.csv"))

# Cdkn1a+ DEG---------
sample_name = sample_ls[1]
modules = module_ls[[sample_name]]
modules = modules[grepl("Quie",names(modules))]
data_S = get(paste0("data_S_",sample_name))
genes = c("Cdkn1a","Lef1")
df = FetchData(data_S, vars = c(genes, modules))
thred = 0.4
df$status = ifelse(data_S@meta.data[,modules] > thred, names(modules), 
                   ifelse(df$Cdkn1a > mean(df$Cdkn1a), "Cdkn1a+", "Cdkn1a-"))
write.csv(df,file = file.path(folder.path,"e13_smom2_cdkn1a_table.csv"))

# Proportion of CC-phase in each module --------
# proportions_df <- df %>%
#   group_by(wnt_related_module) %>%
#   summarise(
#     G1_count = sum(Phase == "G1"),  # Count of G1 cells in each stage
#     total_count = n(),  # Total number of cells in each stage
#     proportion = G1_count / total_count  # Proportion of G1 cells
#   ) %>%
#   ungroup()
# ggplot(proportions_df, 
#        aes(x=wnt_related_module, y=proportion, group = 1)) +
#   geom_line() + 
#   geom_point(aes(size = total_count)) + 
#   labs(size = "Cell number", color = "Condition") + 
#   ggtitle("Proportion of cells in G1 Phase") + 
#   theme_minimal()

