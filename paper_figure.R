# Sensitivity Test----------
## Rank Stability-----------
score_df = lapply(c(3:10,20,30,40,50),function(kNN){
  try({
    W = Symmetric_KNN_graph(knn = kNN, feature_space = feature_space)[[1]]
    rho = Rowwise_normalize(dat[,colnames(W)])
    marker = fast_calculate_multi_score(W = W, init_state = rho)
    
    sub_score0 = grep("score0",colnames(marker))
    sub_score1 = grep("score1_0",colnames(marker))
    accu_score1 = apply(marker,1,function(x){
      sum(x[sub_score0]/x[sub_score1])
    })
    accu_score1
  },silent = TRUE)
})
names(score_df) = paste0("kNN",c(3:10,20,30,40,50))
score_df = do.call(cbind,score_df)
score_ranked <- data.frame(score_df) %>% mutate(across(everything(), rank))
base = "kNN5"
top_genes_list <- lapply(c(50,seq(100,1000,100)), 
                         function(x) head(score_ranked[order(score_ranked[,base]),], x))
average_ranks <- lapply(top_genes_list, function(top_genes) {
  colMeans(top_genes)
})
average_ranks <- data.frame(
  Top_Genes = c(50,seq(100,1000,100)),
  do.call(rbind, average_ranks)
)
average_ranks <- average_ranks %>% tidyr::pivot_longer(cols = -Top_Genes, names_to = "kNN", values_to = "Average_Rank")
average_ranks$kNN = as.numeric(gsub("kNN","",average_ranks$kNN))
ggplot(data = average_ranks, aes(x = kNN, y = Average_Rank, colour = factor(Top_Genes), group = Top_Genes)) +
  geom_line() +
  geom_point() +
  xlab("Top Genes") +
  ylab(sprintf("Average Rank of Top # Genes\n (base on %s)",base)) +
  theme_minimal() + facet_wrap(~Top_Genes,scales = "free")

## AUC stability------------
gt_genes_list <- c(lapply(c(50,seq(100,400,100)), function(x){
  names(max_logfc)[1:x]
}),list(true_marker))
auc_df = apply(score_ranked,2,function(x){
  unlist(lapply(gt_genes_list,function(genes){
    roc = pROC::roc(as.integer(rownames(score_ranked) %in% genes), x, direction = ">")
    pROC::auc(roc)
  }))
})
auc_df = data.frame(auc_df,gt = c(paste0("top",c(50,seq(100,400,100))),"cellMarkerDB"))
auc_df <- auc_df %>% tidyr::pivot_longer(cols = -gt, names_to = "kNN", values_to = "AUC")
auc_df$kNN = as.numeric(gsub("kNN","",auc_df$kNN))
ggplot(data = auc_df, aes(x = kNN, y = AUC, colour = factor(gt,levels = unique(gt)), group = gt)) +
  geom_line() +
  geom_point() +
  xlab("kNN") +
  ylab("AUC") +
  theme_minimal()
