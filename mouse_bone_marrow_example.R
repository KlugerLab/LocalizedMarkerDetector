dir.path <- "/banach1/ruiqi/local_marker/LocalMarkerDetector/"
source(file.path(dir.path,"LMD_function.R"))
# Define package lists
bioc_packages <- c("clusterProfiler", "AnnotationDbi", "ReactomePA", "org.Mm.eg.db","gprofiler2", "msigdbr")
github_packages <- c("satijalab/seurat-wrappers");
names(github_packages) = strsplit(github_packages, "/")[[1]][2]
names(github_packages)[grepl("seurat-wrappers",names(github_packages))] = "SeuratWrappers"
# Install packages
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
sapply(1:length(github_packages), function(i) if (!requireNamespace(names(github_packages)[i], quietly = TRUE)) remotes::install_github(github_packages[i]))
lapply(c(bioc_packages,names(github_packages)), require, character.only = TRUE)

# Download/Load Data ===========
dir.path <- "/banach1/ruiqi/local_marker/LMD_data"
folder.path <- file.path(dir.path,"tabular_muris")
dir.create(folder.path, recursive=T)
tissue_download_link = c(
  "marrow_facs" = "https://figshare.com/ndownloader/files/13092380",
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
    tiss <- RunALRA(tiss, assay = "RNA"); DefaultAssay(tiss) = "RNA"
    saveRDS(tiss,file = file.path(folder.path,file_name))
  }
}

# Prepare input data ===========
tissue_name = "marrow_facs"
tiss <- readRDS(file.path(folder.path,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
n_dim = 20
feature_space = as.matrix(tiss@reductions$pca@cell.embeddings[,1:n_dim])
dat = as.matrix(tiss[[DefaultAssay(tiss)]]@data)
Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_genes = (Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)
selected_genes = names(selected_genes)[selected_genes]
dat = dat[selected_genes,,drop = FALSE]

# RunLMD step by step ---------------
dir.create(file.path(folder.path,"intermediate_data",tissue_name), recursive=T)
## Construct KNN graph -----------
res = Symmetric_KNN_graph(knn = 5, feature_space = feature_space)
A = res$adj_matrix # Adjacency Matrix
W = res$graph # Symmetrized Graph Matrix
rm(res)
## Diffused graph -------
P_ls = Obtain_Pls(W = W, max_time = 2^20)
rho = Rowwise_normalize(dat[,colnames(W)])
res = fast_get_lmds(W = W, init_state = rho, P_ls = P_ls)
saveRDS(P_ls, file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_diffused_graph.rds")))
saveRDS(res, file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_lmd_profile.rds")))
saveRDS(W, file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_knn_graph.rds")))
saveRDS(rho, file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_gene_distribution.rds")))

# Identify Gene Modules---------
local_gene = show_result_lmd(res)
dat_alra = as.matrix(tiss[["alra"]]@data)[names(local_gene$cut_off_gene),]
dist = Calculate_distance(dat_alra, method = "jaccard")
res = Obtain_gene_partition(dist, clustering_method = "average", 
                            deepSplit = 1, return_tree = TRUE)
saveRDS(local_gene,file = file.path(folder.path,"intermediate_data",tissue_name,paste0("local_gene_",tissue_name,".rds")))
saveRDS(dist, file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_gene_dist.rds")))
saveRDS(res, file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_gene_tree.rds")))

## Module Activity Scores ---------
tiss <- readRDS(file.path(folder.path,paste0(tissue_name,".rds")))
local_gene = readRDS(file.path(folder.path,"intermediate_data",tissue_name,paste0("local_gene_",tissue_name,".rds")))
res = readRDS(file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_gene_tree.rds")))
tiss = AddModuleActivityScore(tiss, gene_partition = res$gene_partition, do_local_smooth = FALSE)
# Filter-out modules express in less than 10 cells
cell_block = tiss@meta.data[,grep("Module",colnames(tiss@meta.data))]
tiss@meta.data = tiss@meta.data[,!grepl("Module",colnames(tiss@meta.data))]
modules = names(which(colSums(cell_block > 0.5) >= 10))
modules = modules[order(as.numeric(gsub("Module","",modules)))]
# Rename modules
modules_rename = setNames(paste0("Module",1:length(modules)),modules)
cell_block = cell_block[,names(modules_rename)]
colnames(cell_block) = modules_rename
tiss = AddMetaData(tiss,cell_block)

res$gene_partition_rename = res$gene_partition
levels(res$gene_partition_rename) = gsub("Module","",modules_rename[paste0("Module",levels(res$gene_partition_rename))])
res$gene_partition_rename = res$gene_partition_rename[!is.na(res$gene_partition_rename)]
local_gene$'gene_partition' = res$gene_partition_rename

saveRDS(tiss,file = file.path(folder.path,paste0(tissue_name,".rds")))
saveRDS(local_gene, file.path(folder.path,"intermediate_data",tissue_name,paste0("local_gene_",tissue_name,".rds")))
saveRDS(res, file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_gene_tree.rds")))

## Project CC modules from SmoM2 onto this embedding
local_gene <- readRDS("/banach1/ruiqi/Peggy_data/Peggy_scdata/240329_smom2/intermediate_data/local_gene_smom2_dermal_E13.5_MUT.rds")
gene_partition = droplevels(local_gene$gene_partition[local_gene$gene_partition %in% c(1,2,3)])
levels(gene_partition) = paste0("smom2_",levels(gene_partition))
tiss = AddModuleActivityScore(tiss, gene_partition = gene_partition, do_local_smooth = FALSE)
FeaturePlot(tiss, features = paste0("Module",levels(gene_partition)), order = TRUE, reduction = "tsne", ncol = 3) & NoAxes() & 
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) & labs(color = "ModuleScore") & NoLegend()
FeaturePlot(tiss, features = paste0("Module",levels(gene_partition)), order = TRUE, reduction = "tsne_cc", ncol = 3) & NoAxes() & 
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) & labs(color = "ModuleScore") & NoLegend()

# Module-Celltype-alignment ---------
# Create one-hot matrix of Cell Type
cell_label = as.factor(tiss$cell_ontology_class)
abbrev = levels(cell_label)
abbrev[grepl("T cell",abbrev)] = "T cell" # combine T cells
abbrev[grepl("natural killer",abbrev)] = "NK cell" # combine NK cells
abbrev[abbrev == "B cell"] = "Cd3e+ Klrb1+ B cell"
abbrev[grepl("Slamf1-positive",abbrev)] = "Slamf1+ MPC"
abbrev[grepl("Slamf1-negative",abbrev)] = "Slamf1- MPC"
abbrev[grepl("common lymphoid",abbrev)] = "CLP"
abbrev[grepl("megakaryocyte",abbrev)] = "MEP"
abbrev[grepl("hematopoietic",abbrev)] = "HSC"
abbrev[grepl("hematopoietic",abbrev)] = "HSC"
abbrev[grepl("granulocyte monocyte",abbrev)] = "GMP"
names(abbrev) = levels(cell_label)
levels(cell_label) = abbrev[levels(cell_label)]
one_hot_matrix <- sapply(unique(cell_label), function(x) as.integer(cell_label == x))
colnames(one_hot_matrix) = unique(cell_label)
rownames(one_hot_matrix) = names(cell_label)
saveRDS(one_hot_matrix, file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_ct_onehot_mtx.rds")))

# Calculate Jaccard index between module activity & cell type
jaccard_index <- apply(tiss@meta.data[,grepl("Module",colnames(tiss@meta.data))], 2, function(module_col) {
  apply(one_hot_matrix, 2, function(one_hot_col) {
    1 - philentropy::distance(rbind(module_col,one_hot_col), method = "jaccard", mute.message=TRUE)
  })
})
saveRDS(jaccard_index, file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_ct_module_jaccard.rds")))

# Transfer to Droplet Dataset ----------
## Module Activity Score
tissue_name = "marrow_droplet"
tiss <- readRDS(file.path(folder.path,paste0(tissue_name,".rds")))
DefaultAssay(tiss) <- "RNA"
tiss = AddModuleActivityScore(tiss, gene_partition = res$gene_partition_rename, do_local_smooth = TRUE)
saveRDS(tiss,file = file.path(folder.path,paste0(tissue_name,".rds")))

## Calculate Jaccard index between module activity & cell type
dir.create(file.path(folder.path,"intermediate_data",tissue_name),recursive = TRUE)
cell_label = tiss$cell_ontology_class
cell_label[!cell_label %in% names(abbrev)] = NA
cell_label = setNames(abbrev[cell_label],names(cell_label))
one_hot_matrix <- sapply(unique(na.omit(cell_label)), function(x) as.integer(cell_label == x))
one_hot_matrix[is.na(one_hot_matrix)] = 0
saveRDS(one_hot_matrix, file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_ct_onehot_mtx.rds")))
jaccard_index <- apply(tiss@meta.data[,grepl("Module",colnames(tiss@meta.data))], 2, function(module_col) {
  apply(one_hot_matrix, 2, function(one_hot_col) {
    1 - philentropy::distance(rbind(module_col,one_hot_col), method = "jaccard", mute.message=TRUE)
  })
})
saveRDS(jaccard_index, file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_ct_module_jaccard.rds")))

# Pathway over-representation Analysis ------------
tissue_name = "marrow_facs"
local_gene = readRDS(file.path(folder.path,"intermediate_data",tissue_name,paste0("local_gene_",tissue_name,".rds")))
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
saveRDS(epathway_result,file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_reactome_result.rds")))

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
saveRDS(ego_result,file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_gomf_result.rds")))

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
saveRDS(msig_result,file = file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_msig_result.rds")))

# Module Annotation table--------
tissue_name = "marrow_facs"
local_gene = readRDS(file.path(folder.path,"intermediate_data",tissue_name,paste0("local_gene_",tissue_name,".rds")))
tiss = readRDS(file.path(folder.path,paste0(tissue_name,".rds")))

gene_partition = local_gene$gene_partition
df_module = data.frame(Module = levels(gene_partition))
df_module[,"Gene"] = unlist(lapply(levels(gene_partition),function(x){paste(names(gene_partition)[gene_partition == x],collapse = ",")}))
df_module[,"# of Genes"] = unlist(lapply(levels(gene_partition),function(x){sum(gene_partition == x)}))
df_module[,"# of Expressed Cells"] = colSums(tiss@meta.data[,paste0("Module",levels(gene_partition))] > 0.5)
# CellType-specific modules
jaccard_index = readRDS(file.path(folder.path,"intermediate_data",tissue_name,paste0(tissue_name,"_ct_module_jaccard.rds")))
thred = 0.4
module_celltype = apply(jaccard_index,2,function(x){paste(names(which(x >= thred)),collapse = ",")})
df_module[,"celltype"] = module_celltype[paste0("Module",levels(gene_partition))]
# Annotate each module with Top5 significant Terms (p.adjust < 0.05)
for(term in c("gomf","reactome","msig")){
  e_result = readRDS(file = file.path(folder.path,"intermediate_data",tissue_name,sprintf("%s_%s_result.rds",tissue_name,term)))
  top_5_term = lapply(levels(gene_partition),function(i){
    terms = e_result[[i]]@result
    terms = terms %>% filter(p.adjust < 0.05) %>% 
      arrange(p.adjust) %>% slice_head(n = 5) %>%.$Description
    paste0(terms,collapse = ", ")
  })
  df_module[,paste0("Top5 Significant ",term)] = unlist(top_5_term)
}
write.csv(df_module,file = file.path(folder.path,"intermediate_data",tissue_name, "module_description.csv"),row.names = FALSE)

# Annotate Cell Cycle ------
mmus_s = gprofiler2::gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gprofiler2::gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
for(tissue_name in c("marrow_facs","marrow_droplet")){
  tiss = readRDS(file.path(folder.path,paste0(tissue_name,".rds")))
  tiss = CellCycleScoring(tiss, s.features = mmus_s, 
                          g2m.features = mmus_g2m, set.ident = FALSE)
  tmp <- tiss %>%
    ScaleData(verbose = FALSE, features = c(mmus_s,mmus_g2m)) %>% 
    RunPCA(npcs = 20, verbose = FALSE, features = c(mmus_s,mmus_g2m)) %>% 
    RunTSNE(dims = 1:20, seed.use = 42, reduction.name = "tsne_cc")
  tiss[["tsne_cc"]] = tmp[["tsne_cc"]]
  rm(tmp)
  saveRDS(tiss, file = file.path(folder.path,paste0(tissue_name,".rds")))
}

## Reactome/Msig table --------
term = "reactome"; tissue_name = "marrow_facs"
cc_module = c(22,20,19)
e_result = readRDS(file = file.path(folder.path,"intermediate_data",tissue_name,sprintf("%s_%s_result.rds",tissue_name,term)))
e_result = e_result[cc_module]
# Annotate each modules with top5 Pathway terms
top_enrich_pathway_cc <- do.call(rbind,lapply(1:length(e_result), function(i) data.frame(e_result[[i]]@result,module = names(e_result)[i])))
top_enrich_pathway_cc <- top_enrich_pathway_cc[complete.cases(top_enrich_pathway_cc), ]
top_enrich_pathway_cc = top_enrich_pathway_cc %>% filter(p.adjust < 0.05)
write.csv(top_enrich_pathway_cc,file = file.path(folder.path, "intermediate_data",tissue_name,"reactome_pathway_cc.csv"))

library(ReactomeGraph4R)
login()
matchObject(id = 'R-MMU-1640170')
matchHierarchy(id = "R-MMU-1640170", type = "row")
