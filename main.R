#####                 NETWORK TAB              #####
create_network_graph <- function(mat) {
  g <- graph_from_adjacency_matrix(mat, 
                                   weighted = T, 
                                   mode = "undirected")
  vnames <- colnames(mat) 
  if (substr(colnames(mat)[1], 1, 4) == "ENSG") {
    vnames <- lapply(colnames(mat), function(n) {
      n <- substr(n, 5, nchar(n)) %>% as.numeric() 
    })
  }
  vertex_attr(g) <- list(
    color = rep("seagreen2", gorder(g)), 
    name = colnames(mat),                      
    label = vnames,                     # Vertex labels
    label.cex = rep(1, gorder(g)),      
    label.color = rep("black", gorder(g)), 
    label.family = rep("sans", gorder(g)), 
    label.dist = rep(2, gorder(g))      # Label distance
  )
  
  return(g)
}

# subset counts matrix with input genes
subset_by_genes <- function(data, genes) {
  rownames(data) <- data[,1]
  if (length(genes) == 0) {
    return(data[-1])
  }
  sbst <- data %>% 
    filter(.[[1]] %in% genes)
  rownames(sbst) <- sbst[,1]
  return(sbst[-1])
}


#####                 DE TAB              #####
plot_de_volcano <- function(data, fc_threshold, padj_threshold, c1, c2) {
  # Ensure padj is non-zero
  data <- data %>% mutate(padj = ifelse(padj == 0, 1e-300, padj))
  thresh <- 10 ^ padj_threshold
  # Classify points
  data <- data %>% mutate(volc_plot_status = case_when(
    padj < thresh & log2FoldChange < (fc_threshold * -1) ~ "DOWN",
    padj < thresh & log2FoldChange > fc_threshold ~ "UP",
    .default = "NS"))

  # Set factor levels explicitly, allowing for empty levels
  data$volc_plot_status <- factor(data$volc_plot_status, levels = c("DOWN", 'NS', "UP"))

  g <- ggplot(data) +
  geom_point(aes(x = log2FoldChange,
                 y = -log10(padj),
                 colour = volc_plot_status)) +
  theme_minimal() +
  labs(x = expression(Log[2]("Fold Change")), 
       y = expression(-Log[10](italic("padj"))),
       title = "Volcano Plot of DESeq2 Differential Expression Results") + 
  scale_color_manual(name = "Differential Expression Status",
                     values = c("DOWN" = c2,  "UP" = c1, "NS" = "grey"), 
                     labels = c("Downregulated", "Upregulated", "Not significant"), 
                     drop = F) +
  theme(legend.position = "bottom", 
        plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  
  return(g)
}

plot_de_log2fc <- function(data, padj_threshold = 0.10) {
  g <- data %>%
    filter(padj < padj_threshold) %>% 
    ggplot() +
    geom_histogram(aes(x = log2FoldChange), 
                   fill = "seagreen2", 
                   color = "black", 
                   bins = 100) +
    theme_minimal() + 
    labs(x = expression(Log[2]("Fold Change")), 
         y = "Count",
      title = "Histogram of Log2FoldChange for DE Genes (control vs. Huntington's)") +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  return(g)
}

plot_de_pvals <- function(data) {
  g <- data %>%
    ggplot() +
    geom_histogram(aes(x = pvalue), 
                   color = "black", 
                   fill = "forestgreen",
                   bins = 50) +
    theme_minimal() + 
    labs(x = "Raw p-values", 
         y = "Count",
         title = "Histogram of raw pvalues obtained from DE analysis (control vs. Huntington's)") +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  return(g)
}

#####                 COUNTS TAB              #####
plot_counts_pca <- function(data, first, second) {
  pca_results <- get_pca_results(data[-1], x = T)
  vars <- get_pca_results(data[-1], x = F)
  pca_results["diagnosis"] <- rep(c("C","H"), c(49, 20))
  g <- ggplot(pca_results) +
    geom_point(aes(x = !!sym(first),
                   y = !!sym(second), 
                   color = !!sym("diagnosis"))) +
    labs(x = glue("{first} ({round(vars[[first]], 2)}% variance explained)"),
         y = glue("{second} ({round(vars[[second]], 2)}% variance explained)"),
         title = glue("{first} vs. {second}")) +
    theme_classic() + 
    scale_color_manual(name = "Diagnosis",
                      values = c("C" = "seagreen2", "H" = "forestgreen"), 
                      labels = c("Neurologically Normal", "Huntington's Disease")) + 
    theme(legend.position = "bottom", 
          plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  return(g)
}

get_pca_results <- function(counts, x = T) {
  pca_results <- prcomp(t(counts), center = T, scale = F) 
  pcs <- data.frame(pca_results$x)
  if (!x) {
    vars <- pca_results$sdev ** 2 
    perc_vars <- vars / sum(vars) * 100
    names(perc_vars) <- colnames(pcs)
    return(perc_vars)
  }
  return(pcs)
}

plot_counts_scatter <- function(counts, val) {
  counts$keep <- factor(counts$keep, levels = c("TRUE", "FALSE"))
  g <- ggplot(counts) + 
    theme_classic() +
    scale_color_manual(name = "Gene Passes Filters", 
                       values = c("TRUE" = "forestgreen", 'FALSE' = "seagreen2"),
                       labels = c("TRUE" = "True", "FALSE" = "False"), 
                       drop = F)
  if (val == 1) {
    g <- g + geom_point(aes(x = log2(!!sym("medians") + 1), 
                            y = log10(!!sym("variance") + 1), 
                            color = !!sym("keep"))) +
      labs(x = "Log2(Median)", y = "Log10(Variance)", 
           title = "Median Count vs. Variance")
  } else {
    g <- g + geom_point(aes(x = log2(!!sym("medians") + 1), 
                            y = !!sym("num_zeros"), 
                            color = !!sym("keep"))) +
      labs(x = "Log2(Median)", y = "Number of Zeros", 
           title = "Median Count vs. Number of Zeros")
  }
  g <- g + theme(legend.position = "bottom",
                 plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  return(g)
}

process_counts_filters <- function(counts, min_nonzeros, var_percentile, names = 0) {
  # filter by nonzeros
  filtered <- counts[rowSums(counts[-1] != 0) >= min_nonzeros,] 
  #filter by variance
  filtered["variance"] <- apply(filtered[-1], MARGIN = 1, FUN = var, na.rm = T)
  var_threshold <- var_percentile / 100 * max(filtered$variance, na.rm = T)
  filtered <- filtered[filtered$variance <= var_threshold, ]
  filtered["variance"] <- NULL
  if (names == 0) {
    return(filtered)
  } else {
    return(rownames(filtered))
  }
 
}

process_counts_summary <- function(counts, filtered) {
  counts <- counts[-1] # get rid of GeneID column
  results <- data.frame(matrix(ncol = 2, nrow = 4))
  results[, 1] <- c("Total samples", "Total genes", "Genes passing filter", "Genes not passing filter")
  colnames(results) <- NULL
  tot_genes <- nrow(counts)
  filtered_genes <- nrow(filtered) 
  results[, 2] <- c(ncol(counts), tot_genes, 
                    glue("{filtered_genes} ({round(filtered_genes/tot_genes*100, 2)}%)"), 
                    glue("{tot_genes -  filtered_genes} ({round((tot_genes - filtered_genes) / tot_genes * 100, 2)}%)"))
  return(results)
}

#####                 SAMPLES TAB              #####
plot_samples_scatter <- function(samples, xcol, ycol, xax, yax, title, colore = NULL) {
  g <- ggplot(samples) + theme_classic()
  if (!is.null(colore)) {
    sbst <- samples %>%
    filter(!is.na(!!sym(colore)))
    g <- g + geom_point(data = sbst, 
                        aes(x = !!sym(xcol), y = !!sym(ycol), 
                            color = as_factor(!!sym(colore)))) +
      scale_color_manual(values = c("forestgreen","seagreen2"), 
                         name = fix_name(colore)) + 
      theme(legend.position = "bottom", 
            plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  } else {
    g <- g + geom_point(aes(x = !!sym(xcol), y = !!sym(ycol)), 
                        color = "forestgreen") +
      theme(plot.title = element_text(hjust=0.5, face = "bold", size = 13))
  }
  g <- g + labs(x = xax, y = yax, title = title) 
  return(g)
}

plot_samples_density <- function(samples, xax) {
  g <- ggplot(samples) +
    geom_density(aes(x = !!sym(xax), fill = !!sym("diagnosis")), alpha = 0.6) +  
    scale_fill_manual(name = "Diagnosis",
                       values = c("Neurologically_normal" = "seagreen2", "Huntington's_Disease" = "forestgreen"), 
                       labels = c("Neurologically Normal", "Huntington's Disease")) +
    labs(x = fix_name(xax), 
         y = "Density",
         title = glue("Density Plot of Patients' {fix_name(xax)}")) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 13),
          legend.position = "bottom")
  return(g)
}

plot_samples_boxplot <- function(samples, selected_col) {
  both <- names(samples)[colSums(is.na(samples)) < 49]
  g <- ggplot(samples)
  if (selected_col %in% both) {
    g <- g +
      geom_boxplot(aes(x = diagnosis, y = !!sym(selected_col), fill = diagnosis)) + 
      scale_fill_manual(values = c("seagreen2", "forestgreen")) +
      scale_x_discrete(labels = c("Huntington's_Disease" = "Huntington's Disease", 
                                  "Neurologically_normal" = "Neurologically Normal"))
  } else {
      sbst <- samples %>% 
        filter(!is.na(!!sym(selected_col)))
      g <- ggplot(sbst, aes(x = diagnosis, y=!!sym(selected_col))) +
        geom_violin(width = 0.4, fill = "forestgreen") +
        geom_boxplot(width = 0.1, fill = "seagreen2") + 
        geom_jitter(color= "black", size = 0.8) +
        scale_x_discrete(labels = c("Huntington's_Disease" = "Huntington's Disease"))
  }
  g <- g + guides(fill = "none") + 
    theme_classic() + 
    labs(x = "Diagnosis", y = fix_name(selected_col), 
         title = glue("{fix_name(selected_col)} in Patient Samples")) +
    theme(plot.title = element_text(hjust=0.5, face = "bold", size = 13))

  return(g)
}  

process_samples_summary <- function(samples) {
  # calculate means and sds of numeric values
  means <- samples %>% summarise_if(is.numeric, mean, na.rm = TRUE) %>% round()
  sds <- samples %>% summarise_if(is.numeric, sd, na.rm = TRUE) %>% round(2)
  vals <- rbind(means, sds) %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "Column") %>% 
    mutate(Values = paste0(V1, " (+/- ", V2, ")")) %>%
    select(c(Column, Values))
  
  # combine values into summary table
  summary_table <- as.data.frame(colnames(samples))
  colnames(summary_table)[1] <- "Column"
  summary_table <- summary_table %>% 
    mutate("Type" = sapply(samples, FUN = class)) %>% 
    left_join(vals, join_by(Column))
  summary_table$Values[1] <- paste(levels(samples$tissue), collapse = ", ")
  summary_table$Values[2] <- paste(levels(samples$diagnosis), collapse = ", ")
  summary_table$Column <- fix_name(summary_table$Column, single = F)
  return(summary_table)
}

# makes column / axes / title names capitalized and remove underscores
fix_name <- function(n, single = T) {
  if (single) {
    x <- str_to_title(gsub("_", " ", n))
  } else {
    x <- lapply(n, FUN = function(y) {
      fix_name(y)
    })
  }
  return(x)
}
