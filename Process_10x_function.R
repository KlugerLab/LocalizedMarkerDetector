Load10xData <- function(data_dir, sample_names, sub_rna_dir = "filtered_feature_bc_matrix", with.protein = F, with.tcr = F, sub_tcr_dir = "TCR", human_only = F){
  data_S_list <- list()
  data_dir0 <- data_dir
  if (!with.protein){
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data_S_list[[i]] <- CreateSeuratObject(Read10X(data.dir = data_dir), project = i)
    }
  }
  if (with.protein){
    # For output from CellRanger >= 3.0 with multiple data types
    #list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
    data_dir0 <- data_dir
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_rna_dir, sep = "/")
      print(data_dir)
      data <- Read10X(data.dir = data_dir)
      data_S_list[[i]] <- CreateSeuratObject(counts = data$`Gene Expression`, project = i)
      data_S_list[[i]][["ADT"]] <- CreateAssayObject(counts = data$`Antibody Capture`)
    }
  }
  if (with.tcr){
    for(i in sample_names){
      data_dir <- paste(data_dir0, i, sub_tcr_dir, sep = "/")
      print(data_dir)
      data_S_list[[i]] <- add_clonotype(data_S_list[[i]], data_dir)
    }
  }
  if (human_only){
    for(i in names(data_S_list)){
      data_S_list[[i]] <- data_S_list[[i]][grep("^hg19", rownames(data_S_list[[i]])), ]
    }
  }
  if (TRUE){
    if (length(grep("^mt",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^mt"
    if (length(grep("^MT",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^MT"
    if (length(grep("^hg19-MT",rownames(data_S_list[[1]]))) > 0) pattern_mt <- "^hg19-MT"
    print(pattern_mt)
    for(i in names(data_S_list)){
      data_S_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
        object = data_S_list[[i]], pattern = pattern_mt
      )
    }
  }
  
  data_S_list
}
