add_newline_in_middle <- function(input_string, max_length) {
  if (nchar(input_string) > max_length) {
    # Find the approximate middle
    middle <- ceiling(nchar(input_string) / 2)
    
    # Find the closest space to the middle
    space_before <- regexpr(" ", substring(input_string, 1, middle), fixed = TRUE)
    space_after <- regexpr(" ", substring(input_string, middle + 1), fixed = TRUE)
    
    if (space_before == -1 && space_after == -1) {
      # No spaces found, cannot split without breaking words
      return(input_string)
    }
    
    if (space_before == -1) {
      split_pos <- middle + space_after
    } else if (space_after == -1) {
      split_pos <- space_before
    } else {
      split_pos <- ifelse(middle - space_before < space_after, space_before, middle + space_after)
    }
    
    # Insert newline at the split position
    output_string <- paste0(
      substring(input_string, 1, split_pos), "\n", substring(input_string, split_pos + 1)
    )
  } else {
    output_string <- input_string
  }
  
  return(output_string)
}

Visualize_jaccard_mtx <- function(jaccard_mtx, threshold = 0.4){
  if(is.null(rownames(jaccard_mtx))|is.null(colnames(jaccard_mtx))){
    stop("Please define celltype names & module names")
  }
  df = reshape2::melt(jaccard_mtx)
  colnames(df)[1:2] = c("CellType","Module")
  df$CellType = factor(df$CellType,levels = rownames(jaccard_mtx))
  df$Module = factor(df$Module,levels = rev(colnames(jaccard_mtx)))
  
  p <- ggplot(df, aes(x=CellType, y=Module, fill=value)) +
    geom_tile(color="white") +
    scale_fill_gradientn(colors=colorRampPalette(rev(c("firebrick3", "white")))(99),limits = c(0,1), name="Jaccard Index") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1,size = 15),
      axis.text.y = element_text(size = 15),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 15),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) + geom_text(data = df %>% filter(value >= threshold), 
                  aes(label = sprintf("%.2f", value)), 
                  color = "white", size = 4)
  return(p)
}

findCommonPrefix <- function(strings) {
  prefix <- strings[1]
  for (str in strings) {
    while (substr(str, 1, nchar(prefix)) != prefix) {
      prefix <- substr(prefix, 1, nchar(prefix) - 1)
      if (nchar(prefix) == 0) break
    }
  }
  return(prefix)
}

