########################### was in UI DE TAB ########################### 
uiOutput("de_tab")
# plotOutput("de_jitter"),

########################### was in server DE TAB ########################### 
# renders DE table sorting parameters once data loads
output$de_tab <- renderUI({
  data = de_data()
  req(data)
  cols = colnames(data)
  sidebarLayout(
    sidebarPanel(
      selectInput("de_col_sortby",
                  label = "Sort by:",
                  choices = cols,
                  selected = cols[1]),
      radioButtons("de_asc_desc",
                  label = "Ascending or descending?",
                  choices = c("Asc", "Desc"),
                  selected = "Asc"),
      numericInput("de_nrows",
                   label = "Number of rows to display",
                   min = 1,
                   max = nrow(data),
                   value = 6),
      actionButton("de_table_btn",
                   label = "Sort Table")
    ),
    mainPanel(
      DT::dataTableOutput("mytable")
    )
)})

# reactive variables for updating DE table
de_table_reactives <- reactiveValues(sort = NULL,
                                            asc = "Asc",
                                            row = 6)
# update DE table by selected parameters
observeEvent(input$de_table_btn, {
  de_table_reactives$sort = input$de_col_sortby
  de_table_reactives$asc = input$de_asc_desc
  de_table_reactives$row = input$de_nrows
})

# outputs DE table
output$de_table <- renderTable({
  req(de_data())
  table <- process_table_sorting(de_data(),
                                 de_table_reactives$sort,
                                 de_table_reactives$asc,
                                 de_table_reactives$row)
  return(table)
})

# outputs scatter plot of top 10 most DE genes by padj
# output$de_jitter <- renderPlot({
#   req(de_data())
#   g <- plot_de_jitter(de_data())
#   return(g)
# })

########################### was in UI DE TAB ########################### 
uiOutput("samples_table_layout")

########################### was in UI SAMPLES TAB ########################### 
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

########################### was in MAIN ########################### 
# plot_de_jitter(data) {
#   top10 <- data %>% 
#     arrange(desc(padj)) %>% 
#     slice_head(n = 10)
#   g <- ggplot(top10) + 
#     geom_point(aes(x = "Gene", 
#                    y = ""), 
#                position = "jitter")
#   return(g)
# }



process_table_sorting <- function(data, sort_col, asc, row) {
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