options(shiny.maxRequestSize = 100 * 1024^2)

library(shiny)
library(DT)
library(LPEseq2)

ui <- fluidPage(
  titlePanel("LPEseq2: Local Pooled Error-Based ANOVA"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Upload input files"),
      
      fileInput(
        "counts_file",
        "Upload counts CSV",
        accept = c(".csv")
      ),
      
      fileInput(
        "meta_file",
        "Upload metadata CSV",
        accept = c(".csv")
      ),
      
      tags$hr(),
      
      h4("2. Select analysis options"),
      
      uiOutput("group_var_ui"),
      
      selectInput(
        "normalize_method",
        "Normalization method",
        choices = c("library_size", "TMM", "DESeq2", "none"),
        selected = "library_size"
      ),
      
      checkboxInput(
        "log_transform",
        "Log2 transform",
        value = TRUE
      ),
      
      numericInput(
        "min_count",
        "Minimum count",
        value = 5,
        min = 0
      ),
      
      numericInput(
        "prior_count",
        "Prior count",
        value = 1,
        min = 0
      ),
      
      numericInput(
        "n_bin",
        "Number of bins",
        value = 100,
        min = 5
      ),
      
      numericInput(
        "df",
        "Spline degrees of freedom",
        value = 10,
        min = 2
      ),
      
      selectInput(
        "trim_method",
        "Trimming method",
        choices = c("mad", "quantile", "fixed"),
        selected = "mad"
      ),
      
      numericInput(
        "d",
        "Fixed trimming threshold",
        value = 1.2,
        min = 0.01
      ),
      
      checkboxInput(
        "use_weighted_between",
        "Use weighted between-group differences",
        value = FALSE
      ),
      
      helpText(
        "If checked, between-group differences are also used for variance trend estimation. ",
        "This can be useful when replicate information is limited, but may inflate variance ",
        "when many genes are truly differentially expressed."
      ),
      
      selectInput(
        "p_method",
        "P-value method",
        choices = c("chisq", "F_inf"),
        selected = "chisq"
      ),
      
      tags$hr(),
      
      actionButton(
        "run",
        "Run LPE-ANOVA",
        class = "btn-primary"
      ),
      
      br(),
      br(),
      
      downloadButton(
        "download_results",
        "Download results"
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Instructions",
          h4("Input format"),
          p("Counts file: genes as rows and samples as columns."),
          p("Metadata file: samples as rows and variables as columns."),
          p("The column names of the counts file must match the row names of the metadata file."),
          tags$hr(),
          h4("Example counts CSV"),
          verbatimTextOutput("counts_example"),
          h4("Example metadata CSV"),
          verbatimTextOutput("meta_example")
        ),
        
        tabPanel(
          "Counts preview",
          DTOutput("counts_preview")
        ),
        
        tabPanel(
          "Metadata preview",
          DTOutput("meta_preview")
        ),
        
        tabPanel(
          "Results",
          DTOutput("results_table")
        ),
        
        tabPanel(
          "Log",
          verbatimTextOutput("log_text")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$counts_example <- renderText({
    "gene,sample1,sample2,sample3,sample4
gene1,100,120,80,95
gene2,50,60,55,70
gene3,10,15,30,28"
  })
  
  output$meta_example <- renderText({
    "sample,group
sample1,Control
sample2,Control
sample3,Treatment
sample4,Treatment"
  })
  
  counts_data <- reactive({
    req(input$counts_file)
    
    counts <- read.csv(
      input$counts_file$datapath,
      row.names = 1,
      check.names = FALSE
    )
    
    counts <- as.matrix(counts)
    
    storage.mode(counts) <- "numeric"
    
    counts
  })
  
  meta_data <- reactive({
    req(input$meta_file)
    
    read.csv(
      input$meta_file$datapath,
      row.names = 1,
      check.names = FALSE
    )
  })
  
  output$counts_preview <- renderDT({
    req(counts_data())
    
    datatable(
      head(counts_data(), 20),
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  output$meta_preview <- renderDT({
    req(meta_data())
    
    datatable(
      meta_data(),
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  output$group_var_ui <- renderUI({
    req(meta_data())
    
    selectInput(
      "group_var",
      "Group variable",
      choices = colnames(meta_data()),
      selected = colnames(meta_data())[1]
    )
  })
  
  analysis_result <- eventReactive(input$run, {
    counts <- counts_data()
    meta <- meta_data()
    
    validate(
      need(!is.null(rownames(meta)), "Metadata must have sample names as row names."),
      need(!is.null(colnames(counts)), "Counts must have sample names as column names."),
      need(
        setequal(colnames(counts), rownames(meta)),
        "Sample names do not match between counts columns and metadata row names."
      ),
      need(input$group_var %in% colnames(meta), "Selected group variable is not in metadata.")
    )
    
    design_formula <- as.formula(paste("~", input$group_var))
    
    prep <- LPE_preprocess(
      counts = counts,
      colData = meta,
      design = design_formula,
      normalize.method = input$normalize_method,
      log.transform = input$log_transform,
      min.count = input$min_count,
      prior.count = input$prior_count,
      verbose = FALSE
    )
    
    res <- LPE_ANOVA(
      object = prep,
      n.bin = input$n_bin,
      df = input$df,
      trim.method = input$trim_method,
      d = input$d,
      use_weighted_between = input$use_weighted_between,
      verbose = FALSE,
      p.method = input$p_method
    )
    
    res
  })
  
  output$results_table <- renderDT({
    req(analysis_result())
    
    datatable(
      analysis_result(),
      options = list(
        scrollX = TRUE,
        pageLength = 20
      )
    )
  })
  
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("LPEseq2_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(analysis_result(), file, row.names = FALSE)
    }
  )
  
  output$log_text <- renderPrint({
    cat("LPEseq2 web tool\n")
    cat("1. Upload counts CSV.\n")
    cat("2. Upload metadata CSV.\n")
    cat("3. Select group variable.\n")
    cat("4. Click Run LPE-ANOVA.\n")
    cat("\n")
    cat("Counts columns must match metadata row names.\n")
    cat("use_weighted_between: ", input$use_weighted_between, "\n")
  })
}

shinyApp(ui = ui, server = server)