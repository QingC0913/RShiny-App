process_counts_filters <- function(counts, min_nonzeros, var_percentile) {
  print(head(counts))
  # filter by nonzeros
  filtered <- counts[rowSums(counts != 0) >= min_nonzeros,] 
  #filter by variance
  filtered["variance"] <- apply(filtered, MARGIN = 1, FUN = var, na.rm = T)
  print(min_nonzeros)
  print(var_percentile)
  print(filtered)
  var_threshold <- var_percentile / 100 * max(filtered$variance, na.rm = T)
  filtered <- filtered %>% 
    filter(variance <= var_threshold) 
  
  return(filtered)
}
process_counts_summary <- function(counts, filtered) {
  results <- data.frame(matrix(ncol = 2, nrow = 4))
  results[, 1] <- c("Total samples", "Total genes", "Genes passing filter", "Genes not passing filter")
  colnames(results) <- NULL
  tot_genes <- nrow(counts)
  filtered_genes <- nrow(filtered) 
  results[, 2] <- c(ncol(counts), tot_genes, 
                    glue("{filtered_genes} ({filtered_genes/tot_genes*100}%)"), 
                    glue("{tot_genes -  filtered_genes} ({(tot_genes - filtered_genes) / tot_genes * 100}%)" ))
  return(results)
}
  
  
plot_samples_scatter <- function(samples, xcol, ycol, colore = NULL) {
  g <- ggplot(samples) + theme_linedraw()
  if (!is.null(colore)) {
    sbst <- samples %>%
    filter(!is.na(!!sym(colore)))
    g <- g + geom_point(data = sbst, 
                        aes(x = !!sym(xcol), y = !!sym(ycol), 
                            color = as_factor(!!sym(colore))), 
                            size = 3) +
      scale_color_manual(values = c("#8402ab","darkturquoise"), 
                         name = colore) +
      theme(legend.position = "bottom")
} else {
    g <- g + geom_point(aes(x = !!sym(xcol), y = !!sym(ycol)), 
                        color = "darkturquoise",  size = 3)

  }
  # g <- g 
  return(g)
}

plot_samples_boxplot <- function(samples, selected_col) {
  both <- names(samples)[colSums(is.na(samples)) < 49]
  g <- ggplot(samples)
  if (selected_col %in% both) {
    g <- g +
      geom_boxplot(aes(x = diagnosis, y = !!sym(selected_col), fill = diagnosis)) + 
      scale_fill_manual(values = c("#2C5F2D", "#97BC62")) 
  } else {
      sbst <- samples %>% 
        filter(!is.na(!!sym(selected_col)))
      g <- ggplot(sbst, aes(x = diagnosis, y=!!sym(selected_col))) +
        geom_violin(width = 0.4, fill = "#2C5F2D") +
        geom_boxplot(width=0.1, fill = "#97BC62") + 
        geom_jitter(color= "black", size = 0.8)
  }
  g <- g + guides(fill = "none") + 
    theme_linedraw()

  return(g)
}  

process_samples_table <- function(data, sort_col, asc, row) {
  if (is.null(sort_col)) {
    sort_col <- colnames(data)[1]
  }
  if (asc == "Asc") {
    data <- data %>% arrange(!!sym(sort_col)) 
  }
  else {
    data <- data %>% arrange(desc(!!sym(sort_col))) 
  }
  return(head(data, row))
} 

process_samples_summary <- function(samples) {
  # calculate means and sds of numeric values
  means <- samples %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% round()
  sds <- samples %>% summarise_if(is.numeric, sd, na.rm = TRUE) %>% round(2)
  vals <- rbind(means, sds) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Column") %>% 
    mutate(Val = paste0(V1, " (+/- ", V2, ")")) %>%
    select(c(Column, Val))
  
  # combine values into summary table
  summary_table <- as.data.frame(colnames(samples))
  colnames(summary_table)[1] <- "Column"
  summary_table <- summary_table %>% 
    mutate("Type" = sapply(samples, FUN = class)) %>% 
    left_join(vals, join_by(Column))
  summary_table$Val[1] <- paste(levels(samples$tissue), collapse = ", ")
  summary_table$Val[2] <- paste(levels(samples$diagnosis), collapse = ", ")
  return(summary_table)
}