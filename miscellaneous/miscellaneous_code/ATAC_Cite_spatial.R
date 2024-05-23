library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(cowplot)
library(dplyr)
source("/banach1/ruiqi/local_marker/LocalMarkerDetector/LMD_function.R", echo = F)

# CITE-seq: 25 proteins ----------
InstallData("bmcite")
bm <- LoadData(ds = "bmcite")
DefaultAssay(bm) <- 'ADT'
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
bm <- RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')
DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
        repel = TRUE, label.size = 2.5) + NoLegend()
FeaturePlot(bm, features = c("adt_CD45RA","adt_CD16","adt_CD161"),
            reduction = 'adt.umap', max.cutoff = 2, ncol = 3)

# ATAC-seq ----------------
library(SeuratData)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
# BiocManager::install("Bioconductor/GenomeInfoDb")

## Load Data & QC -------
# InstallData("pbmcMultiome")
pbmc.atac <- LoadData("pbmcMultiome", "pbmc.atac")
# QC
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered" & nCount_ATAC < 7e4 & nCount_ATAC > 5e3)
# Extract_t_cells
annotated_cell_types = names(table(pbmc.atac$seurat_annotations))
sub.pbmc.atac <- subset(pbmc.atac, seurat_annotations %in% annotated_cell_types[!grepl("B|HSPC|Mono|Plasma|DC|NK",annotated_cell_types)])

## Process -------
tiss <- sub.pbmc.atac
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
Annotation(tiss) <- annotations
tiss <- RunTFIDF(tiss)
tiss <- FindTopFeatures(tiss, min.cutoff = "q0")
tiss <- RunSVD(tiss)
tiss <- RunUMAP(tiss, reduction = "lsi", dims = 2:50, 
                     reduction.name = "umap.atac", reduction.key = "atacUMAP_")
DimPlot(tiss, group.by = "seurat_annotations", label = TRUE) + 
  NoLegend() + ggtitle("ATAC")

## LMD -------
feature_space = Embeddings(tiss@reductions$lsi)[,2:50]
dat = as.matrix(tiss[["ATAC"]]@data)
Peak_detected_count <- apply(dat > apply(dat,2,median),1,sum)
selected_peaks = (Peak_detected_count >= 10) & (Peak_detected_count <= ncol(dat) * 0.5)
selected_peaks = names(selected_peaks)[selected_peaks]
dat = dat[selected_peaks,,drop = FALSE]

tic()
lmd_result = LMD(dat, feature_space, max_time = 2^20, highres = FALSE)
toc()
local_gene = show_result_lmd(lmd_result)

FeaturePlot(tiss,names(local_gene$gene_rank)[1:50],order = TRUE, ncol = 10) & NoAxes() & NoLegend()

# Group into Modules --------------
tiss <- RunALRA(tiss)
DefaultAssay(tiss) = "ATAC"
dat_alra = as.matrix(tiss[["alra"]]@data)[names(local_gene$cut_off_gene),]
distance_method = "jaccard"
dist = Calculate_distance(dat_alra, method = distance_method)
clustering_method = "average"
gene_hree = hclust(dist, method = clustering_method)
gene_partition = dynamicTreeCut::cutreeDynamic(
  dendro = gene_hree, 
  distM = as.matrix(dist), deepSplit = 2,
  pamStage = TRUE, pamRespectsDendro = TRUE,
  minClusterSize = 10)
names(gene_partition) = labels(dist)
gene_partition = as.factor(gene_partition)
module_order = order(sapply(levels(gene_partition), function(mod) {
  median(which(gene_partition[gene_hree$order] == mod))
}))
levels(gene_partition)[module_order] = seq(nlevels(gene_partition))
gene_partition = factor(gene_partition,levels = seq(nlevels(gene_partition)))
local_gene$'gene_partition' = gene_partition

# ModuleScore -----------
tiss = AddModuleActivityScore(tiss, local_gene$gene_partition, 
                                do_local_smooth = TRUE, assay = "ATAC")
pl_module = FeaturePlot_meta(tiss,feature_partition = local_gene$gene_partition,
                             reduction = "umap.atac", assays = "ATAC")
(wrap_plots(pl_module,ncol = 10) & 
    labs(color = "AvgExp") & 
    NoAxes())
pl_module = lapply(levels(local_gene$gene_partition),function(module){
  FeaturePlot(tiss, features = paste0("Module",module), order = TRUE) + 
    # scale_color_gradientn(colours = rev(brewer_pal(palette = "RdYlBu")(10)), limits = c(0,1)) +
    scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0,1)) +
    ggtitle(sprintf("Module %s (%d genes)",module, sum(local_gene$gene_partition == module)))
})
names(pl_module) = levels(local_gene$gene_partition)
(wrap_plots(pl_module,ncol = 10) & 
    labs(color = "ModuleScore") & 
    NoAxes()) + plot_layout(guides = "collect")
module_df = data.frame(module = levels(local_gene$gene_partition),
                       min_rank = unlist(lapply(levels(local_gene$gene_partition),function(i) min(local_gene$gene_rank[names(local_gene$gene_partition)[local_gene$gene_partition == i]]))),
                       cell_num = colSums(tiss@meta.data[,paste0("Module",levels(local_gene$gene_partition))] > 0.5),
                       gene_num = as.numeric(table(local_gene$gene_partition)) )
modules = module_df$module[module_df$gene_num >= 10 & 
                             (module_df$min_rank < 500 & 
                                module_df$cell_num > 5)]
ClosestFeature(tiss, regions = names(local_gene$gene_partition)[local_gene$gene_partition == 57])
# Module-wise Top10 LMD genes & DoHeatmaps ------
top10 = do.call(rbind,lapply(modules,function(module){
  genes = names(local_gene$gene_partition)[local_gene$gene_partition == module]
  rank = local_gene$gene_rank[genes]
  genes = genes[order(rank)]
  if(min(rank) >= 500) return(NULL)
  else{return(data.frame("gene" = genes[1:10], "module" = module))}
}))
tiss <- ScaleData(tiss, assay = "alra", features = top10$gene)
DoHeatmap(tiss, features = top10$gene, 
          group.by = "seurat_annotations", 
          slot = "scale.data", assay = "alra",
          label = TRUE) + 
  theme(axis.text.y = element_text(size = 0)) + scale_fill_viridis() 

# Compare with Wilcox.rank.sum ---------
# marker <- FindAllMarkers(tiss, features = selected_peaks, 
#                          only.pos = FALSE, test.use = "wilcox")
marker <- presto:::wilcoxauc.Seurat(X = subset(tiss,features = selected_peaks), 
                                   assay = "data", 
                                   seurat_assay = "ATAC")
colnames(marker)[1] = "gene"

# marker$lmd_rank = local_gene$cut_off_gene[marker$gene]
# marker %>%
#   group_by(cluster) %>%
#   filter(p_val_adj < 0.05) %>%
#   slice_head(n = 10) %>%
#   ungroup() -> top10
df_label = data.frame("lmd" = local_gene$gene_rank[selected_peaks], 
                      "wilcox.no.filter" = match(selected_peaks,marker %>% arrange(.,padj,desc(logFC)) %>% distinct(.,gene, .keep_all = TRUE) %>%.$gene),
                      label = selected_peaks)
df_label[c(which(df_label[,1] <= 1000 & df_label[,2] > 10000)),"color"] = colnames(df_label)[1]
df_label[c(which(df_label[,1] > 10000 & df_label[,2] <= 1000) ),"color"] = colnames(df_label)[2]
df_label$'color' = factor(df_label$'color',levels = colnames(df_label)[1:2])
rank_cor = stats::cor(log(df_label[,1:2],base = 10),method = "pearson")[1,2]
df_label %>% filter(color == "lmd") %>% arrange(.,lmd)
g = c("chr14-105817737-105818927","chr10-48671056-48673240",
      "chr19-50723222-50725949","chr19-50657772-50659599")
ggplot(df_label, aes(x=get(colnames(df_label)[2]), y=get(colnames(df_label)[1]),color = color)) + geom_point(size = 1) + 
  scale_color_manual(values = c("red","blue","grey50"),
                     labels = c("Favored by LMD","Favored by Wilcox")) + 
  geom_abline(intercept = 0, slope = 1, col = "grey", linetype = "dashed") + 
  scale_x_log10(breaks = 10^(0:4), labels = expression(10^0, 10^1, 10^2, 10^3, 10^4)) +
  scale_y_log10(breaks = 10^(0:4), labels = expression(10^0, 10^1, 10^2, 10^3, 10^4)) + 
  expand_limits(x = c(1, 10^4), y = c(1, 10^4)) +
  labs(x = colnames(df_label)[2],
       y = colnames(df_label)[1],
       color = "Highlighted Genes",
       title = paste("cor =", round(rank_cor, 2))) + 
  ggrepel::geom_text_repel(data = df_label %>% filter(label %in% g), aes(label=label),min.segment.length = unit(0, 'lines'),nudge_y = -0.5,show.legend = FALSE,size = 5)

FeaturePlot(tiss,features = g,order = TRUE) & NoLegend() & NoAxes()
# Spatial (Visium, sequencing-based)