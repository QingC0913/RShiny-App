#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
source("main.R")
library(shiny)
library(tidyverse)
library(rlang) # make sure col names can be used in shiny
library(glue) # string concats
library(DT) # datatable
library(igraph) # network analysis 
library(pheatmap) # heatmap
library(stringr) # string capitalization
library(colourpicker) # color inputs


options(shiny.maxRequestSize = 30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  includeCSS("www/style.css"),
    titlePanel("mRNA-Seq Expression Profiling of Human Post-Mortem BA9 Brain Tissue for Huntington's Disease and Neurologically Normal Individuals"), 
     tabsetPanel(
       tabPanel("Sample Information Exploration",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("samples_file",
                              "Please upload a sample information CSV file.",
                              multiple = F,
                              accept = ".csv"),
                    actionButton(inputId = "samples_upload_btn", 
                                 label = span("Upload File", 
                                 icon("open-file", lib = "glyphicon"))), 
                    width = 3),
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Data Information", 
                               htmlOutput("samples_info")),
                      tabPanel("Sample Summary",
                               textOutput("samples_summary"),
                               tableOutput("samples_table_summary")),
                      tabPanel("Data Table",
                               dataTableOutput("samples_table_layout")
                               ),
                      tabPanel("Sample Exploration Plots",
                               uiOutput("samples_plots_layout")))))),
       tabPanel("Counts Matrix Exploration",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("counts_file",
                              label = "Please upload a normalized RNA-seq counts CSV file.",
                              multiple = F,
                              accept = ".csv"),
                    actionButton("counts_upload_btn", 
                                 label = span("Upload File", 
                                 icon("open-file", lib = "glyphicon"))),
                    uiOutput("counts_param"), 
                    width = 3),
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Data Summary",
                               tableOutput("counts_summary")),
                      tabPanel("Diagnostic Plots",
                               plotOutput("counts_var_plot1"),
                               plotOutput("counts_var_plot2")),
                      tabPanel("Counts Heatmap",
                               plotOutput("counts_heatmap")),
                      tabPanel("Principal Component Analysis",
                               uiOutput("counts_pca_layout")
                               ))))),
       tabPanel("Differential Expression Analysis",
                sidebarLayout(
                  sidebarPanel(
                    fileInput("de_file",
                              label = "Please upload a differential expression results CSV file.",
                              accept = ".csv",
                              multiple = F),
                    actionButton("de_btn",
                                 label = span("Upload File", 
                                 icon("open-file", lib = "glyphicon"))), 
                    uiOutput("de_ui"),
                    width = 3),
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Differential Expression Results",
                               dataTableOutput("de_table")
                               ),
                      tabPanel("Visualizations", 
                               plotOutput("de_volcano"),
                               plotOutput("de_pval_hist"), 
                               plotOutput("de_log2fc_hist")))))),
       tabPanel("Correlation Network Analysis", 
        sidebarLayout(
          sidebarPanel(
              fileInput("network_file", 
                        label = "Please upload a normalized RNA-seq counts CSV file.", 
                        multiple = F, 
                        accept = ".csv"),
              actionButton("network_btn", 
                           label = span("Upload File", 
                           icon("open-file", lib = "glyphicon"))), 
              uiOutput("network_ctrls"), 
              width = 3), 
          mainPanel(
                tabsetPanel(
                  tabPanel("Selected Genes Heatmaps", 
                           verbatimTextOutput("genes_not_found"),
                           plotOutput("network_cor_heatmap"),
                           plotOutput("network_heatmap")
                           ), 
                  tabPanel("Correlation Network", 
                           uiOutput("corr_ui")
                           ),
                  tabPanel("Network Metrics", 
                           dataTableOutput("network_metrics"))))))))

server <- function(input, output, session) {
  #####               NETWORK TAB            #####
  
  # gets normalized counts data for network functionalities
  network_data <- eventReactive(input$network_btn, {
    file = input$network_file
    netw <- read.csv(file$datapath)
    return(netw)
  })
  
  # outputs correlation heatmap
  output$network_cor_heatmap <- renderPlot({
    mat <- correlation_mat(rplc = F)
    if (nrow(mat) > 35) {
      h <-  pheatmap(mat, show_rownames = F, show_colnames = F,
                     main = "Heatmap of Pairwise Gene Expression Correlation")
    } else {
      h <- pheatmap(mat, main = "Heatmap of Pairwise Gene Expression Correlation")
    }
    return(h)
  })
  
  # outputs correlation / counts heatmaps
  output$network_heatmap <- renderPlot({
    sbst <- subset_by_genes(network_data(), get_genes())
    if (nrow(sbst) > 35) {
      h <- pheatmap(log10(sbst + 1), show_rownames = F, 
                    main = "Heatmap of Expression of Selected Genes (log10-transformed)")
    } else {
      h <- pheatmap(log10(sbst + 1), 
                    main = "Heatmap of Expression of Selected Genes (log10-transformed)")
    }
    return(h)
  })  
  
  # correlation parameters for network graph
  output$network_ctrls <- renderUI({
    req(network_data())
    tagList(
      br(), 
      textAreaInput("network_genes", 
                    label = "Please enter a list of genes, one gene per line:",
                    placeholder = "all genes"), 
      actionButton("network_genes_btn", 
                   label = span("Upload Gene List", 
                   icon("upload", lib = "glyphicon"))))
  })
  
  # get list of genes from text area
  get_genes <- eventReactive(input$network_genes_btn, {
    genes <- input$network_genes %>% 
      strsplit(split = "\n", fixed = T)
    genes <- genes[[1]] %>% 
                  lapply(str_trim) %>%
                  unlist()
    genes <- as.character(genes[genes != ""])
    return(genes)
  })
  
  # text output if any input genes are not found in expr matrix
  output$genes_not_found <- renderText({
    data <- network_data()
    genes <- get_genes()
    if (length(genes) == 0) {
      return("")
    }
    not_found <- genes[! genes %in% data[,1]]
    if (length(not_found) == 0) {
      return("")
    }
    return(paste0("The following genes were not found in the counts matrix: ", 
                  toString(not_found)))
  })
  
  # outputs controls for network graph node selection
  output$corr_ui <- renderUI({
    genes <- rownames(subset_by_genes(network_data(), get_genes())) %>% sort()
    req(genes)
    sidebarLayout(
    mainPanel(
      plotOutput("network_graph"),
      verbatimTextOutput("network_shortest")
    ),
    sidebarPanel(
      selectInput("net_layout", 
                  label = "Graph Layout Style", 
                  choices = c("Grid", "Tree", "Circle", "Star", "Random", 
                              "Fruchterman-Reingold", "Kamada-Kawai"),
                  selected = "Grid"), 
      sliderInput("network_slider", 
                  label = "Minimum correlation threshold:", 
                  value = 0.5, 
                  min = 0, 
                  max = 1, 
                  step = 0.05), 
        selectInput("node1", 
                    label = "Gene 1", 
                    choices = genes, 
                    selected = genes[1]), 
        selectInput("node2", 
                    label = "Gene 2",
                    choices = genes, 
                    selected = genes[2]), 
        actionButton("sp_btn", 
                     label = span("Find Shortest Path")),
        width = 4), 
    position = "right")
  })

  # correlation matrix 
  correlation_mat <- function(rplc = T) {
    data <- network_data()
    genes <- get_genes()
    sbst <- subset_by_genes(data, genes)
    zero_sd_rows <- apply(sbst, 1, sd) == 0
    sbst <- sbst[!zero_sd_rows, ]
    if (!rplc) {
      mat <- cor(t(sbst), method = "pearson")
    } else {
      # positive values-only correlation matrix, for finding shortest path
      mat <- cor(t(sbst), method = "pearson") %>% abs()
      min_r <- input$network_slider %>% abs()
      mat[mat < min_r] <- 0 # no edge if doesn't meet criteria
    }
    return(mat)
  }
 
  # outputs correlation network graph
  output$network_graph <- renderPlot({
    mat <- correlation_mat()
    g <- create_network_graph(mat)
    layout_style <- switch(input$net_layout, 
                           "Grid" = layout_on_grid(g),
                           "Tree" = layout_as_tree(g), 
                           "Circle" = layout_in_circle(g),
                           "Star" = layout_as_star(g),
                           "Random" = layout_randomly(g), 
                           "Fruchterman-Reingold"= layout_with_fr(g),
                           "Kamada-Kawai" = layout_with_kk(g))
    plot(simplify(g),
         layout = layout_style,
         vertex.size = degree(g, mode="all"),
         asp = 0)
  })
  
  # outputs shortest path between two chosen nodes
  output$network_shortest <- renderText({
    return(find_shortest_path())
  })
  
  # finds shortest path once button is clicked
  find_shortest_path <- eventReactive(input$sp_btn, {
    mat <- correlation_mat() 
    g <- create_network_graph(mat)
    sp <- shortest_paths(g, from = input$node1, 
                         to = input$node2)
    sp <- sp$vpath[[1]]
    
    if (length(sp) == 0) {
      txt <- glue("There is no path from {input$node1} to {input$node2}.")
    } else {
      txt <- paste0("Shortest path between ", input$node1, " and ", input$node2,
                    "\n[start] ", toString(names(sp)), " [end]")
    }
    return(txt)
  })
  
  # outputs table of network metrics
  output$network_metrics <- renderDataTable({
    mat <- correlation_mat() 
    req(mat)
    g <- create_network_graph(mat)
    
    degs <- degree(g)
    close_centrality <- closeness(g)
    btw_centrality <- betweenness(g)
    df <- data.frame(Gene = colnames(mat), 
                        Degree = as.integer(degs), 
                        Closeness = round(close_centrality, 3), 
                        Betweenness = as.integer(btw_centrality))
    return(df)
    
  }, rownames = F)
  #####                 DE TAB              #####
  
  # gets DE data from file
  de_data <- eventReactive(input$de_btn, {
    file <- input$de_file
    de = read.csv(file$datapath)
    colnames(de) <- c("Gene", colnames(de)[2:ncol(de)])
    return(de)
  })
  
  # reactive variables for DE table / volcano plot filtering
  de_reactives <- reactiveValues(fc_thresh = 1.2, 
                                 p_thresh = 0, 
                                 up = "forestgreen", 
                                 down = "seagreen2")
  
  observeEvent(input$de_volcano_btn, {
    de_reactives$fc_thresh <- as.numeric(input$log2_slider)
    de_reactives$p_thresh <- as.numeric(input$padj_slider)
    de_reactives$up <- input$de_upreg_color 
    de_reactives$down <- input$de_downreg_color
  })
  
  # outputs DE results data
  output$de_table <- renderDataTable({
    data <- de_data()
    req(data)
    sbst <- data[data$padj < (10^de_reactives$p_thresh), ]
    sbst <- sbst[abs(sbst$log2FoldChange) > de_reactives$fc_thresh, ]
    return(sbst)},
    options = list(scrollX = T),
    rownames = F)
  
  # renders UI for volcano plot and DE table (filter by padj and log2FC)
  output$de_ui <- renderUI({
    req(de_data())
    tagList(
      br(), 
      colourInput(inputId = "de_upreg_color", 
                  label = "Color for upregulated genes", 
                  value = "forestgreen"),
      colourInput(inputId = "de_downreg_color", 
                  label = "Color for downregulated genes:", 
                  value = "seagreen2"),
      sliderInput(inputId = "log2_slider", 
                  label = "Filter by log2(Fold Change):", 
                  min = 1, 
                  max = 5, 
                  step = 0.1,
                  value = 1.2),
      sliderInput(inputId = "padj_slider", 
                  label = "Filter by magnitude of p-adjusted:", 
                  min = -40, 
                  max = 0, 
                  step = 1,
                  value = 0), 
      actionButton(inputId = "de_volcano_btn", 
                   label = span("Update Values", 
                   icon("pencil", lib = "glyphicon"))))
  })
  
  # outputs raw pvalues histogram
  output$de_pval_hist <- renderPlot({
    req(de_data())
    g <- plot_de_pvals(de_data())
    return(g)
  })
  
  # outputs log2FC histogram, after filtering with padj threshold
  output$de_log2fc_hist <- renderPlot({ 
    req(de_data())
    g <- plot_de_log2fc(de_data()) 
    return(g)
  })
  
  # outputs DE volcano plot
  output$de_volcano <- renderPlot({
    req(de_data())
    g <- plot_de_volcano(de_data(), de_reactives$fc_thresh, de_reactives$p_thresh, 
                         de_reactives$up, de_reactives$down) 
    return(g)
  })
    
  #####               COUNTS TAB            #####
  # get counts data from file
  counts_data <- eventReactive(input$counts_upload_btn, {
    file = input$counts_file
    counts = read.csv(file$datapath)
    genes <- counts$GeneID
    rownames(counts) <- genes
    return(counts)
  })

  # render counts data params when counts data is available
  output$counts_param <- renderUI({
    req(counts_data())
    tagList(
      br(),
      sliderInput("counts_var_slider",
                  label = "Include genes with variance in the top X percentile:",
                  min = 0,
                  value = 100,
                  max = 100),
      sliderInput("counts_nonzero_slider",
                  label = "Include genes expressed in at least X samples:",
                  min = 0,
                  value = 0,
                  max = ncol(counts_data()) - 1, # first col is GeneID
                  step = 1),
      actionButton("counts_btn", 
                   label = span("Set Filters", 
                   icon("filter", lib = "glyphicon"))))
  })
  
  # reactive values to update counts filters
  counts_summ_reactives <- reactiveValues(var_perc = 100,
                                          nonzeros = 0)
  # update counts plots by selected parameters
  observeEvent(input$counts_btn, {
    counts_summ_reactives$var_perc <- input$counts_var_slider
    counts_summ_reactives$nonzeros <- input$counts_nonzero_slider
  })

  # outputs counts summary
  output$counts_summary <- renderTable({
    req(counts_data())
    counts <- counts_data()
    filtered <- process_counts_filters(counts, counts_summ_reactives$nonzeros,
                                       counts_summ_reactives$var_perc, 0)
    results <- process_counts_summary(counts, filtered)
    return(results)
  })

  # outputs counts plot of median counts vs. variance
  output$counts_var_plot1 <- renderPlot({
    req(counts_data())
    counts <- counts_data()
    keep <- process_counts_filters(counts, counts_summ_reactives$nonzeros,
                                   counts_summ_reactives$var_perc, 1)

    v <- apply(counts[-1], MARGIN = 1, FUN = var, na.rm = T)
    m <- apply(counts[-1], MARGIN = 1, FUN = median, na.rm = T)
    k <- rownames(counts) %in% keep

    df <- data.frame(variance = v, medians = m, keep = k)
    g <- plot_counts_scatter(df, 1)
    return(g)
  })

  # outputs counts plot of median counts vs. num of zeros
  output$counts_var_plot2 <- renderPlot({
    req(counts_data())
    counts <- counts_data()
    keep <- process_counts_filters(counts, counts_summ_reactives$nonzeros,
                                   counts_summ_reactives$var_perc, 1)

    n_z <- apply(counts[-1], MARGIN = 1, FUN = function(x) {sum(x == 0)})
    m <- apply(counts[-1], MARGIN = 1, FUN = median, na.rm = T)
    k <- rownames(counts) %in% keep

    df <- data.frame(num_zeros = n_z, medians = m, keep = k)
    g <- plot_counts_scatter(df, 2)
    return(g)
  })

  # reactive values to update PCA plot
  counts_pca_reactives <- reactiveValues(first = "PC1",
                                          second = "PC2")
  # update counts PCA by selected parameters
  observeEvent(input$pca_btn, {
    counts_pca_reactives$first <- input$counts_pca_select1
    counts_pca_reactives$second <- input$counts_pca_select2
  })

  # outputs PCA parameters when PCA results are available
  output$counts_pca_layout <- renderUI({
    req(counts_data())
    pca_results <- get_pca_results(counts_data()[-1], x = T)
    pcs <- colnames(pca_results)
    sidebarLayout(
      mainPanel(
        plotOutput("counts_pca")),
      sidebarPanel(
        selectInput("counts_pca_select1",
                    label = "First PC",
                    choices = pcs,
                    selected = "PC1",
                    multiple = F),
        selectInput("counts_pca_select2",
                    label = "Second PC",
                    choices = pcs,
                    selected = "PC2",
                    multiple = F),
         actionButton("pca_btn",
                      label = span("Plot PCA", 
                      icon("pencil", lib = "glyphicon")))))
  })

  # outputs PCA plot
  output$counts_pca <- renderPlot({
    req(counts_data())
    g <- plot_counts_pca(counts_data(), 
                         counts_pca_reactives$first, 
                         counts_pca_reactives$second)
    return(g)
  })

  # outputs counts heatmap
  output$counts_heatmap <- renderPlot({
    req(counts_data())
    filtered <- process_counts_filters(counts_data(), counts_summ_reactives$nonzeros,
                                       counts_summ_reactives$var_perc, 0)
    mat <- log10(filtered[-1] + 1) %>% as.matrix()
    pheatmap(mat,
             clustering_distance_rows = "euclidean", # Clustering metric for rows
             clustering_distance_cols = "euclidean", # Clustering metric for columns
             clustering_method = "complete",    # Hierarchical clustering method
             show_rownames = F, 
             main = "Heatmap of Normalized & Log10 Transformed RNA-seq Counts")
  })

  #####               SAMPLES TAB            #####
  # outputs information about study
  output$samples_info <- renderUI({
    HTML("
        <div></div>
        <div>
           This dataset is an mRNA expression profile 
           conducted on post-mortem human dorsolateral 
           prefrontal cortex tissue obtained from patients who had 
           Huntington's Disease, alongside age- and sex-matched 
           neurologically healthy control subjects. 
         </div>
         <div></div>
         <div>
           The dataset is available on NCBI's Gene Expression Omnibus 
           (GEO) with the accession code 
           <a href = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810'>
            GSE64810
           </a>.  
         </div>
         <div>
            mRNA expression counts files, including raw sequencing counts, DESeq2-normalized counts, 
          and TPM/FPKM-normalized counts, are available, along with 
          differential expression results analyzed with DESeq2. 
         </div>
         <div></div>
         <div>
            Relevant information to know: 
            <ul>
              <li>
                <strong>PMI</strong>: post-mortem index, refers time interval 
                between patient's time of death and the time that 
                the deceased patient was examined. 
              </li>
              <li>
                 <strong>CAG</strong>: number of CAG-trinucletide repeats in the Huntingtin (HTT) gene.
              </li>
              <li>
                <strong>RIN</strong>: RNA Integrity Number, refers to quality
                of RNA obtained before sequencing was conducted. 
              </li>
              <li>
                <strong>Vonsattel grade</strong>: a measure of severity of 
                neostriatal atrophy, named after Jean-Paul Vonsattel, M.D.
                There are five grades, from 0 (indiscernible neuropathological abnormalities)
                to 4 (severe neuropathological abnormalities).</li>
            </ul>
")
    })
  
  # loads file data when submit button is pressed
  samples_data <- eventReactive(input$samples_upload_btn, {
    file = input$samples_file
    csv = read.csv(file$datapath)
    csv <- csv %>%
      mutate(tissue = as_factor(tissue), diagnosis = as_factor(diagnosis),
             pmi = as.integer(pmi), age_of_death = as.integer(age_of_death),
             rin = as.numeric(rin), mrna_seq_reads = as.integer(mrna_seq_reads),
             age_of_onset = as.integer(age_of_onset), duration = as.integer(duration),
             cag = as.integer(cag), vonsattel_grade = as.integer(vonsattel_grade),
             h_v_striatal_score = as.numeric(h_v_striatal_score),
             h_v_cortical_score = as.numeric(h_v_cortical_score))
    return(csv)
  })

  # outputs samples summary table
  output$samples_table_summary <- renderTable({
    req(samples_data())
    summary_table <- samples_data() %>%
      process_samples_summary() # in main.R
    return(summary_table)
  })

  # outputs samples summary
  output$samples_summary <- renderText({
    req(samples_data())
    return(paste("The samples data contains", nrow(samples_data()), "and",
                 ncol(samples_data()), "columns. Here is a summary:"))
  })

  # reactive values to update samples table
  samples_table_reactives <- reactiveValues(sort = NULL,
                                            asc = "Asc",
                                            row = 6)

  # update samples table by selected parameters
  observeEvent(input$samples_update_table_btn, {
    samples_table_reactives$sort = input$samples_table_radio
    samples_table_reactives$asc = input$samples_table_asc
    samples_table_reactives$row = input$samples_table_rows
  })

  # outputs samples table 
  output$samples_table_layout <- renderDataTable({
    req(samples_data())
    samples <- samples_data() 
    colnames(samples) <- fix_name(colnames(samples), single = F)
    return(samples)},
    options = list(scrollX = T),
    rownames = F)
  
  # reactive var for plotting samples boxplot
  samples_boxplot <- reactiveVal("age_of_death")
  observeEvent(input$samples_boxplot_btn, {
    samples_boxplot(input$samples_box_radio)
  })

  # outputs samples boxplot selectable columns
  output$samples_plots_layout <- renderUI({
    req(samples_data())

    samples  <- samples_data()
    col_names <- names(samples)[sapply(samples, is.numeric)]
    names(col_names) <- fix_name(col_names)
    col_names <- sort(col_names)
    sidebarLayout(
      mainPanel(plotOutput("samples_boxplot"),
                plotOutput("samples_density"),
                plotOutput("samples_point1"),
                plotOutput("samples_point2"),
                plotOutput("samples_point3")),
      sidebarPanel(
        selectInput("samples_box_radio",
                     label = "Choose a column to plot:",
                     choices = col_names,
                     selected = "age_of_death"),
        actionButton("samples_boxplot_btn",
                     label = span("Plot Boxplot", 
                     icon("pencil", lib = "glyphicon")))))
  })

  output$samples_density <- renderPlot({
    req(samples_data())
    g <- plot_samples_density(samples_data(), "age_of_death")
    return(g)
  })
  
  output$samples_point3 <- renderPlot({
    req(samples_data())
    g <- plot_samples_density(samples_data(), "pmi")
    return(g)
  })
  
  # outputs samples boxplot
  output$samples_boxplot <- renderPlot({
    req(samples_data())
    samples <- samples_data()
    g <- plot_samples_boxplot(samples, samples_boxplot())
    return(g)
  })

  output$samples_point1 <- renderPlot({
    req(samples_data())
    g <- plot_samples_scatter(samples_data(), 
                              xcol = "cag", ycol = "age_of_death", 
                              xax = "Number of CAG Tri-nucleotide Repeats", yax = "Age of Death",
                              title = "Number of CAG Repeats vs. Age of Death in Patient Samples")
    return(g)
  })
  output$samples_point2 <- renderPlot({
    req(samples_data())
    g <- plot_samples_scatter(samples_data(), 
                              xcol = "cag", ycol = "age_of_onset", 
                              xax = "Number of CAG Tri-nucleotide Repeats", yax = "Age of Onset",
                              title = "Number of CAG Repeats vs. Age of Onset in Patient Samples")
    return(g)
  })
  
    output$samples_point3 <- renderPlot({
      req(samples_data())
      g <- plot_samples_scatter(samples_data(), xcol = "h_v_striatal_score",
                                ycol = "h_v_cortical_score", colore = "vonsattel_grade", 
                                xax = "Striatal Score", yax = "Cortical Score", 
                                title = "Striatal vs. Cortical Score in Patient Samples")
      return(g)
  })

} # end server


# Run the application
shinyApp(ui = ui, server = server)