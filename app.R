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
                      tabPanel("tab1", 
                               # tableOutput("placeholder"),
                               tableOutput("counts_summary")), 
                      tabPanel("tab2"), 
                      tabPanel("tab3"), 
                      tabPanel("tab4")
                    )
                  )
                )), 
       tabPanel("DE"), 
       tabPanel("TODO")
      )
)

server <- function(input, output) {
  #####               COUNTS TAB            #####
  counts_data <- eventReactive(input$counts_upload_btn, {
    file = input$counts_file
    counts = read.csv(file$datapath) 
    return(counts[1:10, 1:11])
  })
  
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
  
  counts_summ_reactives <- reactiveValues(var_perc = NULL,    
                                          nonzeros = NULL)
  observeEvent(input$counts_btn, {
    counts_summ_reactives$var_perc <- input$counts_var_slider
    counts_summ_reactives$nonzeros <- input$counts_nonzero_slider
  })
  
  output$counts_summary <- renderTable({
    req(counts_data())
    
    counts <- counts_data()
    
    min_nonzeros <- counts_summ_reactives$nonzeros 
    if (is.null(min_nonzeros)) {
      return(counts)
    }
    # exclude first column which is geneID
    filtered <- counts[rowSums(counts[-1] != 0) >= min_nonzeros,] 

    return(filtered)
  })
  
  output$placeholder <- renderTable({
    req(counts_data())
    return(counts_data())
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
