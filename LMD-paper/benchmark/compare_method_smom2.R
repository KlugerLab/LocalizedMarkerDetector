cran_packages <- c("readxl","dplyr","BiocManager","RColorBrewer")
bioc_packages <- c("Seurat")
sapply(cran_packages, function(pkg) if(!requireNamespace(pkg, quietly = TRUE)){install.packages(pkg)})
sapply(bioc_packages, function(pkg) if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg))
lapply(c(cran_packages,bioc_packages), require, character.only = TRUE)
source(file.path(dirname(rstudioapi::getActiveDocumentContext()$path),"run_methods_function.R"))

# Download/Load Data ===========
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../paths.R")

# Load data ===========
dir.path <- dir.path0
data_source = "smom2"
folder.path <- file.path(dir.path,data_source)
sample_ls = unlist(lapply(c("E13.5", "E14.5"),function(i) paste0(i,"_", c("MUT", "CTL"))))
sample_ls = paste0("smom2_dermal_",sample_ls)

data_S_ls = lapply(sample_ls, function(sample_name){
  tiss <- readRDS(file.path(folder.path,"process_data",paste0("data_S_",sample_name,".rds")))
  tiss
})
names(data_S_ls) = sample_ls

# Run ===========
method_ls = c("lmd","hvg", "wilcox_no_filter",
              "haystack","hotspot",
              "cosg","scmarker","semi")

lapply(sample_ls, function(sample_name){
  # Prepare input data --------
  tiss = data_S_ls[[sample_name]]
  DefaultAssay(tiss) <- "RNA"
  cat(sample_name,":",FindPC(tiss),"PC\n")
  dat = as.matrix(tiss[["RNA"]]@data)
  Gene_detected_count <- apply(dat > apply(dat,2,median),1,sum)
  selected_genes = names(Gene_detected_count)[(Gene_detected_count >= 10) & (Gene_detected_count <= ncol(dat) * 0.5)]
  dat = dat[selected_genes,,drop = FALSE]
  feature_space = Embeddings(tiss[["pca"]])[,1:FindPC(srat = tiss)]
  
  # Run Each Method and save results ========
  dir.path <- dir.path0
  folder.path = file.path(dir.path,"benchmark",data_source)
  if(!file.exists(folder.path)){
    dir.create(folder.path, recursive = T)
  }
  
  for(method in method_ls){
    dir.file = file.path(folder.path, paste0(method,"_",sample_name,".csv"))
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
  }
  
  # Load Rank Table ==========
  dir.path <- dir.path0
  folder.path.rank <- file.path(dir.path,"benchmark",data_source)
  
  df_benchmark = lapply(method_ls, function(method){
    marker = read.table(file.path(folder.path.rank, paste0(method,"_",sample_name,".csv")),header = TRUE)
    if(!"gene" %in% colnames(marker)){
      marker$'gene' = rownames(marker)
    }
    marker %>% select(gene,rank) %>% distinct()
  })
  names(df_benchmark) = method_ls
  df_benchmark = bind_rows(df_benchmark, .id = "method")
  df_benchmark <- tidyr::pivot_wider(df_benchmark, names_from = method, values_from = rank) %>% as.data.frame()
  df_benchmark[is.na(df_benchmark)] = nrow(df_benchmark)
  write.table(df_benchmark,file = file.path(folder.path.rank, paste0(sample_name,"_benchmark_rank_table.csv")),row.names = FALSE)
  
})


# Loop over LMD gene modules -------
gene_module_ls = lapply(sample_ls, function(sample_name){
  local_gene = readRDS(file.path(folder.path,"intermediate_data",sprintf("local_gene_%s.rds",sample_name)))
  local_gene$'gene_partition'
})
names(gene_module_ls) = sample_ls
sample_name = sample_ls[1]
df_benchmark = read.table(file.path(dir.path0,"benchmark",data_source, paste0(sample_name,"_benchmark_rank_table.csv")),
                          header = TRUE)
gene_module_ls[sample_name]
pl = lapply(levels(gene_module_ls[[sample_name]]),function(i){
  module = gene_module_ls[[sample_name]]
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
p = wrap_plots(pl,nrow = 4)
ggsave("1.png",p,height = 12, width = 20)

