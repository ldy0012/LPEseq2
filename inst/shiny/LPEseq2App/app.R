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
      
      selectInput(
        "analysis_method",
        "Analysis method",
        choices = c(
          "LPE-ANOVA" = "LPE",
          "Standard one-way ANOVA" = "standard_anova",
          "Auto by group sample size" = "auto"
        ),
        selected = "auto"
      ),
        
      conditionalPanel(
        condition = "input.analysis_method == 'auto'",
        
        numericInput(
          "standard_min_group_n",
          "Minimum group size for standard ANOVA in auto mode",
          value = 5,
          min = 2
        ),
        
        helpText(
          "In auto mode, standard one-way ANOVA is used when every group has at least this number of samples. ",
          "Otherwise, LPE-ANOVA is used."
        )
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
        "Pseudo count",
        value = 1,
        min = 0
      ),
      
      conditionalPanel(
        condition = "input.analysis_method != 'standard_anova'",
        
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
          "Between-group outlier trimming method",
          choices = c(
            "IQR / boxplot rule" = "iqr",
            "None" = "none"
          ),
          selected = "iqr"
        ),
        
        helpText(
          "The IQR method removes boxplot-style outliers from raw between-group log2 differences. ",
          "No user-defined threshold or k value is required. ",
          "Trimming is applied only when between-group differences are used."
        ),
        
        checkboxInput(
          "use_weighted_between",
          "Use weighted between-group differences",
          value = FALSE
        ),
        
        helpText(
          "If checked, between-group differences are also used for variance trend estimation. ",
          "Outlier trimming is applied only to between-group-derived raw log2 differences. ",
          "Within-group differences are retained for variance trend estimation."
        ),
        
        selectInput(
          "p_method",
          "LPE p-value method",
          choices = c("chisq", "F_inf"),
          selected = "chisq"
        )
      ),

      tags$hr(),
      
      actionButton(
        "run",
        "Run Analysis",
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
          "Method info",
          verbatimTextOutput("method_info")
        ),
        
        tabPanel(
          "Trimming info",
          verbatimTextOutput("trim_info"),
          DTOutput("trim_table")
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
    
    validate(
      need(all(is.finite(counts)), "Counts file contains non-numeric, NA, NaN, or Inf values."),
      need(all(counts >= 0), "Counts file contains negative values.")
    )
    
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
    
    meta <- meta[colnames(counts), , drop = FALSE]
    
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
      use_weighted_between = input$use_weighted_between,
      analysis.method = input$analysis_method,
      standard.min.group.n = input$standard_min_group_n,
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
  
  output$method_info <- renderPrint({
    req(analysis_result())
    
    method <- attr(analysis_result(), "analysis.method")
    requested_method <- attr(analysis_result(), "requested.analysis.method")
    standard_min_group_n <- attr(analysis_result(), "standard.min.group.n")
    
    if (is.null(method)) {
      if ("method" %in% colnames(analysis_result())) {
        method <- unique(analysis_result()$method)
      } else {
        method <- "unknown"
      }
    }
    
    if (is.null(requested_method)) {
      requested_method <- input$analysis_method
    }
    
    if (is.null(standard_min_group_n)) {
      standard_min_group_n <- input$standard_min_group_n
    }
    
    cat("Requested analysis method:", requested_method, "\n")
    cat("Actually selected analysis method:", method, "\n")
    cat("Minimum group size for standard ANOVA in auto mode:", standard_min_group_n, "\n")
    
    if (input$analysis_method == "auto") {
      cat("\nAuto mode rule:\n")
      cat("- If every group has at least standard.min.group.n samples: standard one-way ANOVA\n")
      cat("- Otherwise: LPE-ANOVA\n")
    }
  })
  
  output$trim_info <- renderPrint({
    req(analysis_result())
    
    info <- attr(analysis_result(), "trim.info")
    
    if (is.null(info)) {
      method <- attr(analysis_result(), "analysis.method")
      
      cat("No trimming information available.\n")
      
      if (!is.null(method) && method == "standard_anova") {
        cat("This is expected because standard one-way ANOVA does not use LPE variance-trend trimming.\n")
      } else {
        cat("This may occur when no between-group-derived trimming information was produced.\n")
      }
      
      return()
    }
    
    cat("outlier trimming method: ", input$trim_method, "\n")
    cat("Between values before trimming:", info$n_between_before, "\n")
    cat("Between values after trimming:", info$n_between_after, "\n")
    cat("Between values removed:", info$n_between_removed, "\n")
  })
  
  output$trim_table <- renderDT({
    req(analysis_result())
    
    info <- attr(analysis_result(), "trim.info")
    
    if (is.null(info) ||
        is.null(info$threshold.table) ||
        nrow(info$threshold.table) == 0) {
      return(
        DT::datatable(
          data.frame(Message = "No threshold table available.")
        )
      )
    }
    
    DT::datatable(
      info$threshold.table,
      options = list(scrollX = TRUE, pageLength = 10)
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
    cat("4. Click Run Analysis.\n")
    cat("\n")
    cat("Counts columns must match metadata row names.\n")
    cat("analysis method: ", input$analysis_method, "\n")
    
    if (input$analysis_method == "auto") {
      cat("standard.min.group.n: ", input$standard_min_group_n, "\n")
    }
    
    if (input$analysis_method != "standard_anova") {
      cat("use_weighted_between: ", input$use_weighted_between, "\n")
      cat("trimming method: ", input$trim_method, "\n")
      cat("Trimming method:", info$method, "\n")
      cat("Trimming rule:", info$rule, "\n")
      cat("Between values before trimming:", info$n_between_before, "\n")
      cat("Between values after trimming:", info$n_between_after, "\n")
      cat("Between values removed:", info$n_between_removed, "\n")
    }
  })
}

shinyApp(ui = ui, server = server)