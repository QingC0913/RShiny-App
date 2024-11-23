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
library(bslib) # to validate csv files
library(glue) # string concats
library(fields) # heatmap legend
# library(ggbeeswarm) # beeswarm plot

# Define UI for application that draws a histogram
ui <- fluidPage(
     tabsetPanel(
       tabPanel("Samples", 
                sidebarLayout(
                  sidebarPanel(
                    fileInput("samples_file", 
                              "Please upload a RNA sequencing data file.",
                              multiple = F, 
                              accept = ".csv"), 
                    actionButton("samples_upload_btn", "Upload")
                    ), 
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Summary", 
                               textOutput("samples_summary"),
                               tableOutput("samples_table_summary")), 
                      tabPanel("Table", 
                               uiOutput("samples_table_layout")
                               ), 
                      tabPanel("Plots", 
                               uiOutput("samples_plots_layout"))
                    )
                  )
                )
              ), 
       tabPanel("Counts", 
                sidebarLayout(
                  sidebarPanel(
                    fileInput("counts_file", 
                              label = "Please upload a counts matrix file.", 
                              multiple = F, 
                              accept = ".csv"),
                    actionButton("counts_upload_btn", "Upload"), 
                    uiOutput("counts_param")
                  ), 
                  mainPanel(
                    tabsetPanel(
                      tabPanel("Summary", 
                               tableOutput("counts_summary")), 
                      tabPanel("Diagnostic Plots", 
                               # tableOutput("counts_var_plot"),
                               plotOutput("counts_var_plot1"),
                               plotOutput("counts_var_plot2")),
                      tabPanel("Heatmap", 
                               plotOutput("counts_heatmap")), 
                      tabPanel("PCA",
                               uiOutput("counts_pca_layout")
                               )
                      )
                  )
                )), 
       tabPanel("DE"), 
       tabPanel("TODO")
      )
)

server <- function(input, output) {
  #####               COUNTS TAB            #####
  # get counts data from file
  counts_data <- eventReactive(input$counts_upload_btn, {
    file = input$counts_file
    counts = read.csv(file$datapath) 
    genes <- counts$GeneID
    rownames(counts) <- genes
    counts <- counts[1:10, 1:11]
    return(counts)
  })
  
  # render counts data params when counts data is available
  output$counts_param <- renderUI({
    req(counts_data())
    tagList(
      sliderInput("counts_var_slider", 
                  label = "min percentile of variance", 
                  min = 0,
                  value = 100, 
                  max = 100),
      sliderInput("counts_nonzero_slider", 
                  label = "min non-zero samples", 
                  min = 0, 
                  value = 0, 
                  max = ncol(counts_data()) - 1, # first col is GeneID
                  step = 1), 
      actionButton("counts_btn", "Submit")
    )
  })
  
  # reactive values to update counts filters
  counts_summ_reactives <- reactiveValues(var_perc = 100,    
                                          nonzeros = 0)
  observeEvent(input$counts_btn, {
    counts_summ_reactives$var_perc <- input$counts_var_slider
    counts_summ_reactives$nonzeros <- input$counts_nonzero_slider
  })

  # outputs counts summary
  output$counts_summary <- renderTable({
    req(counts_data())
    counts <- counts_data()
    min_nonzeros <- counts_summ_reactives$nonzeros
    var_percentile <- counts_summ_reactives$var_perc
    filtered <- process_counts_filters(counts, min_nonzeros, var_percentile, 0)
    results <- process_counts_summary(counts, filtered)
    return(results)
  })
  
  # outputs counts plot of median counts vs. variance
  output$counts_var_plot1 <- renderPlot({
    req(counts_data())
    counts <- counts_data()
    min_nonzeros <- counts_summ_reactives$nonzeros
    var_percentile <- counts_summ_reactives$var_perc
    keep <- process_counts_filters(counts, min_nonzeros, var_percentile, 1)

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
    min_nonzeros <- counts_summ_reactives$nonzeros
    var_percentile <- counts_summ_reactives$var_perc
    keep <- process_counts_filters(counts, min_nonzeros, var_percentile, 1)

    n_z <- apply(counts[-1], MARGIN = 1, FUN = function(x) {sum(x == 0)})
    m <- apply(counts[-1], MARGIN = 1, FUN = median, na.rm = T)
    k <- rownames(counts) %in% keep

    df <- data.frame(num_zeros = n_z, medians = m, keep = k)
    g <- plot_counts_scatter(df, 2)
    return(g)
  })
  
  output$counts_heatmap <- renderPlot({
    req(counts_data())
    counts <- counts_data()
    min_nonzeros <- counts_summ_reactives$nonzeros
    var_percentile <- counts_summ_reactives$var_perc
    filtered <- process_counts_filters(counts, min_nonzeros, var_percentile, 0)
    mat <- log2(filtered[-1] + 1) %>% as.matrix()
    colors <- colorRampPalette(c("red", "white", "black"))(15)
    
    heatmap(mat, col = colors)
    image.plot(legend.only = T, 
               zlim = range(mat, na.rm = T), 
               col = colors,
               legend.lab = "Normalized Counts")
  })
  
  # todo remove 
  output$counts_pca2 <- renderTable({
    req(counts_data())
    counts <- counts_data()[-1]
    pca_results <- prcomp(t(counts), center = T, scale = F)
    print(pca_results$x)
    
    return(pca_results$x)
  })
  
  counts_pca_reactives <- reactiveValues(first = "PC1",
                                          second = "PC2")
  observeEvent(input$pca_btn, {
    counts_pca_reactives$first <- input$counts_pca_select1
    counts_pca_reactives$second <- input$counts_pca_select2
  })
  
  output$counts_pca_layout <- renderUI({
    req(counts_data())
    pca_results <- get_pca_results(counts_data()[-1], x = T)
    pcs <- colnames(pca_results)
    sidebarLayout(
      sidebarPanel(
        selectInput("counts_pca_select1", 
                    label = "first PC", 
                    choices = pcs, 
                    selected = "PC1", 
                    multiple = F), 
        selectInput("counts_pca_select2", 
                    label = "second PC", 
                    choices = pcs, 
                    selected = "PC2", 
                    multiple = F),
         actionButton("pca_btn", 
                      label = "Plot!")
      ), 
      mainPanel(
        plotOutput("counts_pca")
      ))
  })
  
  output$counts_pca <- renderPlot({
    req(counts_data())
    pca_results <- get_pca_results(counts_data()[-1], x = T)
    vars <- get_pca_results(counts_data()[-1], x = F)
    first <- counts_pca_reactives$first 
    second <- counts_pca_reactives$second
    print(vars)
    g <- ggplot(pca_results) +
      geom_point(aes(x = !!sym(first), 
                     y = !!sym(second))) + 
      labs(x = glue("{first} ({round(vars[[first]], 2)}% variance explained)"), 
           y = glue("{second} ({round(vars[[second]], 2)}% variance explained)"), 
           title = glue("{first} vs. {second}"))
    return(g)
  })
  
  #####               SAMPLES TAB            #####
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
  
  # outputs samples table column radio buttons
  output$samples_table_layout <- renderUI({
    req(samples_data())
    samples <- samples_data()
    col_names = colnames(samples)
    sidebarLayout(
      sidebarPanel(
        radioButtons("samples_table_radio", 
                     label = "Sort by a column", 
                     choices = col_names, 
                     selected = col_names[1]),
        radioButtons("samples_table_asc", 
                     label = "Ascending or descending", 
                     choices = c("Asc", "Desc"), 
                     selected = "Asc"),
        numericInput("samples_table_rows", 
                     label = "Number of rows", 
                     value = 6, 
                     min = 1,
                     max = nrow(samples)),
        actionButton("samples_update_table_btn", "Update Table")), 
      mainPanel(
        tableOutput("samples_table"))
    )})
  
  # outputs samples table
  output$samples_table <- renderTable({
    req(samples_data())
    
    # error if not CSV
    ext <- tools::file_ext(input$samples_file$datapath)
    validate(need(ext == "csv", "Please upload a CSV file"))
    
    # sort / arrange table by parameters
    table <- process_samples_table(samples_data(), 
                                   samples_table_reactives$sort, 
                                   samples_table_reactives$asc, 
                                   samples_table_reactives$row)
    return(table)
  })
  
  # reactive var for plotting samples boxplot
  samples_boxplot <- reactiveVal("age_of_death")
  observeEvent(input$samples_boxplot_btn, {
    samples_boxplot(input$samples_box_radio)
  })
  
  # outputs samples boxplot selectable columns
  output$samples_plots_layout <- renderUI({
    req(samples_data())

    samples  <- samples_data()
    col_names = names(samples)[sapply(samples, is.numeric)]
    sidebarLayout(
      mainPanel(plotOutput("samples_boxplot"), 
                plotOutput("samples_point1"), 
                plotOutput("samples_point2"), 
                plotOutput("samples_point3")),
      sidebarPanel(
        radioButtons("samples_box_radio",
                     label = "Choose a column visualize with boxplot!",
                     choices = col_names,
                     selected = "age_of_death"),
        actionButton("samples_boxplot_btn",
                     label = "Plot boxplot!")
    )
)})

  # outputs samples boxplot
  output$samples_boxplot <- renderPlot({
    req(samples_data())
    samples <- samples_data()
    g <- plot_samples_boxplot(samples, samples_boxplot())
    return(g)
  })
  
  output$samples_point1 <- renderPlot({
    req(samples_data())
    g <- plot_samples_scatter(samples_data(), "cag", "age_of_death")
    return(g)
})
  
  output$samples_point2 <- renderPlot({
    req(samples_data())
    g <- plot_samples_scatter(samples_data(), "age_of_onset", "age_of_death")
    return(g)
})
  
  output$samples_point3 <- renderPlot({
    req(samples_data())
    g <- plot_samples_scatter(samples_data(), "h_v_striatal_score", 
                              "h_v_cortical_score", "vonsattel_grade")
    return(g)
})
  
} # end server
  

# Run the application 
shinyApp(ui = ui, server = server)

#beeswarm 
# pcas <- 9
# pca_long <- pca_results$x[, 1:pcas] %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = paste0("PC", seq(1:pcas)), names_to = "PC", values_to = "vals") 
# 
# print(pca_long)
# g <- pca_long %>% 
#   ggplot + 
#   geom_beeswarm(aes(x = !!sym("PC"), y = !!sym("vals")))