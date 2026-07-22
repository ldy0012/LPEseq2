options(shiny.maxRequestSize = 100 * 1024^2)

library(shiny)
library(DT)
library(LPEseq2)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .run-status-running {
        color: #2563EB;
        font-weight: bold;
        padding: 8px 0;
      }
      .run-status-done {
        color: #16A34A;
        font-weight: bold;
        padding: 8px 0;
      }
      .run-status-error {
        color: #DC2626;
        font-weight: bold;
        padding: 8px 0;
      }
      @keyframes blink {
        0%  { opacity: 1; }
        50% { opacity: 0.3; }
        100%{ opacity: 1; }
      }
      .blinking {
        animation: blink 1s infinite;
      }
    "))
  ),
  titlePanel("LPEseq2: Local Pooled Error-Based ANOVA"),

  sidebarLayout(
    sidebarPanel(
      h4("1. Upload input files"),

      fileInput("counts_file", "Upload counts file", accept = c(".csv", ".tsv", ".txt")),
      checkboxInput(
        "no_gene_id",
        "First column is NOT a gene identifier (auto-assign gene IDs)",
        value = FALSE
      ),
      helpText(
        "Check this if your file has no gene identifier column. ",
        "Gene identifiers will be automatically assigned as gene_1, gene_2, ..."
      ),

      uiOutput("gene_id_warning_ui"),

      fileInput("meta_file", "Upload metadata file", accept = c(".csv", ".tsv", ".txt")),

      tags$hr(),

      h4("2. Select analysis options"),

      uiOutput("group_var_ui"),

      selectInput(
        "normalize_method",
        "Normalization method",
        choices = c("library_size", "TMM", "DESeq2", "none"),
        selected = "TMM"
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
          "Pairwise outlier trimming method",
          choices = c(
            "Pooled bin-wise IQR" = "iqr",
            "Fixed D-value threshold (LPEseq1)" = "dvalue",
            "None" = "none"
          ),
          selected = "iqr"
        ),

        helpText(
          "The IQR method pools within-group and between-group-derived pairwise values, ",
          "divides them into expression-intensity A-bins, and applies the conventional 1.5 within each bin. ",
          "Outlier detection is performed on the M-value scale used for variance trend estimation."
        ),

        conditionalPanel(
          condition = "input.trim_method == 'dvalue'",

          numericInput(
            "d_threshold",
            "D-value threshold",
            value = 1.2,
            min = 0,
            step = 0.1
          ),

          helpText(
            "Applies a single fixed threshold to the raw pairwise difference D, ",
            "as in LPEseq1's non-replicate outlier procedure (Gim et al. 2016). ",
            "Any pairwise value with |D| greater than this threshold is removed, ",
            "regardless of expression-intensity bin. Default 1.2 was empirically ",
            "tuned on specific benchmark datasets; consider adjusting for your data."
          )
        ),

        checkboxInput(
          "use_weighted_between",
          "Use weighted between-group differences",
          value = TRUE
        ),

        helpText(
          "If checked, between-group-derived values are also included in variance trend estimation. ",
          "When IQR trimming is selected, within-group and between-group-derived values are pooled before bin-wise IQR trimming."
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

      conditionalPanel(
        condition = "input.run > 0",
        uiOutput("run_status")
      ),

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
          p("If the first column is not a gene identifier, check 'First column is NOT a gene identifier' to assign gene IDs automatically."),
          tags$hr(),
          h4("Counts file format"),
          verbatimTextOutput("counts_example"),
          h4("Metadata file format"),
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
          "Variance trend info",
          verbatimTextOutput("trend_info"),
          DTOutput("base_var_table")
        ),

        tabPanel(
          "Trimming info",
          verbatimTextOutput("trim_info"),
          DTOutput("trim_table")
        ),

        tabPanel(
          "Spline Plot",
          plotOutput("spline_plot", height = "500px"),
          helpText("Blue dots: bin-level variance estimates | Red line: fitted variance trend spline")
        ),

        # ###
        # tabPanel(
        #   "Volcano Plot",
        #   fluidRow(
        #     column(3,
        #            numericInput("volcano_fc_cutoff", "Mean difference cutoff", value = 1, min = 0, step = 0.1),
        #            numericInput("volcano_q_cutoff", "q-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
        #            helpText("X-axis: Between-group MS (log2 scale) | Y-axis: -log10(p.value)")
        #     ),
        #     column(9,
        #            plotOutput("volcano_plot", height = "500px")
        #     )
        #   )
        # ),
        # ###

        tabPanel(
          "Log",
          verbatimTextOutput("log_text"),
        )
      )
    )
  )
)

server <- function(input, output, session) {

  run_state <- reactiveVal("idle")  # idle / running / done / error
  gene_id_warning <- reactiveVal(NULL)

  observeEvent(input$run, {
    run_state("running")
  })

  output$run_status <- renderUI({
    state <- run_state()
    if (state == "idle") {
      return(NULL)
    } else if (state == "running") {
      div(
        class = "run-status-running blinking",
        icon("spinner"), " Running analysis... Please wait."
      )
    } else if (state == "done") {
      div(
        class = "run-status-done",
        icon("check-circle"), " Analysis complete."
      )
    } else if (state == "error") {
      div(
        class = "run-status-error",
        icon("exclamation-circle"), " An error occurred."
      )
    }
  })

  output$gene_id_warning_ui <- renderUI({
    msg <- gene_id_warning()
    if (is.null(msg)) return(NULL)
    div(
      style = paste(
        "background-color: #FEF9C3;",
        "border-left: 4px solid #EAB308;",
        "padding: 8px 12px;",
        "margin-top: 4px;",
        "font-size: 0.9em;"
      ),
      icon("triangle-exclamation"),
      strong(" Note: "),
      msg
    )
  })

  analysis_result <- eventReactive(input$run, {

    result <- tryCatch({

      counts <- counts_data()
      meta   <- meta_data()

      validate(
        need(!is.null(rownames(meta)), "Metadata must have sample names as row names."),
        need(!is.null(colnames(counts)), "Counts must have sample names as column names."),
        need(!anyDuplicated(colnames(counts)), "Counts sample names must be unique."),
        need(!anyDuplicated(rownames(meta)), "Metadata sample names must be unique."),
        need(
          setequal(colnames(counts), rownames(meta)),
          "Sample names do not match between counts columns and metadata row names."
        ),
        need(input$group_var %in% colnames(meta), "Selected group variable is not in metadata.")
      )

      meta <- meta[colnames(counts), , drop = FALSE]

      validate(
        need(!anyNA(meta[[input$group_var]]), "Selected group variable contains NA values."),
        need(
          length(unique(meta[[input$group_var]])) >= 2,
          "At least two groups are required for analysis."
        )
      )

      design_formula <- stats::reformulate(input$group_var)

      prep <- tryCatch(
        LPE_preprocess(
          counts           = counts,
          colData          = meta,
          design           = design_formula,
          normalize.method = input$normalize_method,
          log.transform    = input$log_transform,
          min.count        = input$min_count,
          prior.count      = input$prior_count,
          verbose          = FALSE
        ),
        error = function(e) {
          stop(paste("Preprocessing failed:", conditionMessage(e)))
        }
      )

      lpe_n_bin              <- if (is.null(input$n_bin)) 100 else input$n_bin
      lpe_df                 <- if (is.null(input$df)) 10 else input$df
      lpe_trim_method        <- if (is.null(input$trim_method)) "iqr" else input$trim_method
      lpe_use_weighted_between <- if (is.null(input$use_weighted_between)) FALSE else input$use_weighted_between
      lpe_d_threshold        <- if (is.null(input$d_threshold)) 1.2 else input$d_threshold
      lpe_p_method           <- if (is.null(input$p_method)) "chisq" else input$p_method
      auto_min_group_n       <- if (is.null(input$standard_min_group_n)) 5 else input$standard_min_group_n

      res <- LPE_ANOVA(
        object             = prep,
        n.bin              = lpe_n_bin,
        df                 = lpe_df,
        trim.method        = lpe_trim_method,
        use_weighted_between = lpe_use_weighted_between,
        d.threshold        = lpe_d_threshold,
        analysis.method    = input$analysis_method,
        standard.min.group.n = auto_min_group_n,
        verbose            = FALSE,
        p.method           = lpe_p_method
      )

      run_state("done")
      res
    }, error = function(e) {
      if (inherits(e, "shiny.silent.error")) {
        run_state("idle")
      } else {
        run_state("error")
      }
      stop(e)
    })
    result
  })

  output$counts_example <- renderText({
    "=== Supported formats ===
- CSV  : comma-separated (.csv)
- TSV  : tab-separated (.tsv, .txt)
- Other: semicolon (;), pipe (|), or space-separated (.txt)
  (separator is detected automatically)

=== Required structure ===
- First column : gene identifiers (row names)
- Other columns: one column per sample (integer counts recommended)
- First row     : header with sample names
- No missing values allowed

=== Example (CSV) ===
gene,sample1,sample2,sample3,sample4
gene1,100,120,80,95
gene2,50,60,55,70
gene3,10,15,30,28

=== Example (TSV) ===
gene    sample1    sample2    sample3    sample4
gene1   100        120        80         95
gene2   50         60         55         70
gene3   10         15         30         28"
  })

  output$meta_example <- renderText({
    "=== Supported formats ===
- CSV  : comma-separated (.csv)
- TSV  : tab-separated (.tsv, .txt)
- Other: semicolon (;), pipe (|), or space-separated (.txt)
  (separator is detected automatically)

=== Required structure ===
- First column : sample names (must match column names of counts file)
- Other columns: one column per variable (e.g. group, batch)
- First row     : header with variable names

=== Example (CSV) ===
sample,group
sample1,Control
sample2,Control
sample3,Treatment
sample4,Treatment

=== Example (TSV) ===
sample    group
sample1   Control
sample2   Control
sample3   Treatment
sample4   Treatment"
  })

  counts_data <- reactive({
    req(input$counts_file)

    counts_raw <- data.table::fread(
      input$counts_file$datapath,
      data.table = FALSE,
      check.names = FALSE
    )

    if (isTRUE(input$no_gene_id)) {
      gene_id_warning(
        paste0(
          "Gene identifiers were automatically assigned as gene_1, gene_2, ... ",
          "(", nrow(counts_raw), " genes total)"
        )
      )
      counts <- as.matrix(counts_raw)
      rownames(counts) <- paste0("gene_", seq_len(nrow(counts_raw)))
    } else {
      gene_id_warning(NULL)
      rownames(counts_raw) <- as.character(counts_raw[[1]])
      counts <- as.matrix(counts_raw[, -1, drop = FALSE])
    }

    storage.mode(counts) <- "numeric"

    validate(
      need(all(is.finite(counts)), "Counts file contains non-numeric, NA, NaN, or Inf values."),
      need(all(counts >= 0), "Counts file contains negative values.")
    )

    counts
  })

  meta_data <- reactive({
    req(input$meta_file)

    meta <- data.table::fread(
      input$meta_file$datapath,
      data.table = FALSE,
      check.names = FALSE
    )
    rownames(meta) <- as.character(meta[[1]])
    meta <- meta[, -1, drop = FALSE]
    meta
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
    trend_info <- attr(analysis_result(), "trend.info")

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

  output$trend_info <- renderPrint({
    req(analysis_result())

    info <- attr(analysis_result(), "trend.info")

    if (is.null(info)) {
      cat("No variance trend information available.\n")
      cat("This is expected when standard one-way ANOVA is selected.\n")
      return()
    }

    cat("Variance trend method:", info$method, "\n")

    if (!is.null(info$spline.df)) {
      cat("Smoothing spline df:", info$spline.df, "\n")
    }
  })

  output$base_var_table <- renderDT({
    req(analysis_result())

    base_var <- attr(analysis_result(), "base.var")

    if (is.null(base_var) || nrow(base_var) == 0) {
      return(
        DT::datatable(
          data.frame(Message = "No bin-level variance points available.")
        )
      )
    }

    DT::datatable(
      base_var,
      options = list(scrollX = TRUE, pageLength = 10)
    )
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
        cat("This may occur when no LPE variance-trend trimming information was produced.\n")
      }

      return()
    }

    cat("Trimming method:", info$method, "\n")
    cat("Trimming rule:", info$rule, "\n")

    if (!is.null(info$trim.scale)) {
      cat("Trimming scale:", info$trim.scale, "\n")
    }

    if (!is.null(info$n_total_before)) {
      cat("Total values before trimming:", info$n_total_before, "\n")
      cat("Total values after trimming:", info$n_total_after, "\n")
      cat("Total values removed:", info$n_total_removed, "\n")
    }

    cat("\nWithin-group values:\n")
    cat("Before:", info$n_within_before, "\n")
    cat("After:", info$n_within_after, "\n")
    cat("Removed:", info$n_within_removed, "\n")

    cat("\nBetween-group-derived values:\n")
    cat("Before:", info$n_between_before, "\n")
    cat("After:", info$n_between_after, "\n")
    cat("Removed:", info$n_between_removed, "\n")
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

  output$spline_plot <- renderPlot({
    req(analysis_result())

    base_var   <- attr(analysis_result(), "base.var")
    trend_info <- attr(analysis_result(), "trend.info")
    var_spline <- attr(analysis_result(), "var.spline")

    if (is.null(base_var) || nrow(base_var) == 0) {
      plot.new()
      text(0.5, 0.5,
           "Spline plot is only available for LPE-ANOVA.\nStandard ANOVA does not estimate a variance trend.",
           cex = 1.2, col = "gray40")
      return()
    }

    x_seq  <- NULL
    y_pred <- NULL

    if (!is.null(var_spline) && var_spline$type == "smooth.spline") {
      x_seq <- seq(min(base_var$A), max(base_var$A), length.out = 300)
      x_seq_clipped <- pmin(pmax(x_seq, var_spline$x_min), var_spline$x_max)
      y_pred <- stats::predict(var_spline$object, x_seq_clipped)$y
    }

    plot(
      base_var$A, base_var$var.M,
      pch  = 16, col = "#3B82F6AA", cex = 0.9,
      xlab = "Mean Expression (A)",
      ylab = "Estimated Local Pooled Variance",
      main = paste0("Intensity-Dependent Variance Trend\n(method: ",
                    if (!is.null(trend_info$method)) trend_info$method else "unknown", ")"),
      las  = 1
    )

    if (!is.null(x_seq) && !is.null(y_pred)) {
      valid <- is.finite(x_seq) & is.finite(y_pred)
      if (any(valid)) {
        lines(x_seq[valid], y_pred[valid], col = "#EF4444", lwd = 2.5)
      }
    }

    legend("topright",
           legend = c("Bin-level variance", "Fitted variance trend"),
           col    = c("#3B82F6AA", "#EF4444"),
           pch    = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2.5),
           bty    = "n")
  })

  # ###
  # output$volcano_plot <- renderPlot({
  #   req(analysis_result())
  #
  #   res <- analysis_result()
  #
  #   # Guard: required columns must exist
  #   if (!all(c("MS_between", "p.value", "q.value") %in% colnames(res))) {
  #     plot.new()
  #     text(0.5, 0.5, "Required columns not found in result.", cex = 1.2, col = "gray40")
  #     return()
  #   }
  #
  #   fc_cut <- input$volcano_fc_cutoff
  #   q_cut  <- input$volcano_q_cutoff
  #
  #   # X-axis: log2-transformed between-group mean square
  #   # Y-axis: -log10 p-value
  #   x <- log2(res$MS_between + .Machine$double.xmin)
  #   y <- -log10(res$p.value  + .Machine$double.xmin)
  #
  #   # Classify genes as significant or not
  #   sig     <- res$q.value < q_cut & res$MS_between > fc_cut
  #   col_vec <- ifelse(sig, "#EF4444", "#94A3B8")
  #   cex_vec <- ifelse(sig, 0.9, 0.7)
  #
  #   plot(
  #     x, y,
  #     col  = col_vec,
  #     pch  = 16,
  #     cex  = cex_vec,
  #     xlab = "log2(MS_between)",
  #     ylab = "-log10(p.value)",
  #     main = paste0("Volcano Plot  (q < ", q_cut, ",  MS_between > ", fc_cut, ")"),
  #     las  = 1
  #   )
  #
  #   # Reference lines for the chosen cutoffs
  #   abline(h = -log10(q_cut),                              col = "#64748B", lty = 2, lwd = 1.2)
  #   abline(v = log2(fc_cut + .Machine$double.xmin),        col = "#64748B", lty = 2, lwd = 1.2)
  #
  #   n_sig <- sum(sig, na.rm = TRUE)
  #   legend("topright",
  #          legend = c(paste0("Significant (n = ", n_sig, ")"),
  #                     paste0("Not significant (n = ", nrow(res) - n_sig, ")")),
  #          col    = c("#EF4444", "#94A3B8"),
  #          pch    = 16, bty = "n", pt.cex = 1)
  #
  #   # Label the top 10 most significant genes
  #   if (n_sig > 0) {
  #     top_idx <- order(res$p.value)[seq_len(min(10, n_sig))]
  #     top_sig <- top_idx[sig[top_idx]]
  #     if (length(top_sig) > 0) {
  #       text(x[top_sig], y[top_sig],
  #            labels = res$gene[top_sig],
  #            cex = 0.65, pos = 3, col = "#1E293B")
  #     }
  #   }
  # })
  # ###

  output$log_text <- renderPrint({
    cat("LPEseq2 web tool\n")
    cat("1. Upload counts file.\n")
    cat("2. Upload metadata file.\n")
    cat("   Note: Counts columns must match metadata row names.\n")  # ← 여기로
    cat("3. Select group variable.\n")
    cat("4. Click Run Analysis.\n")
    cat("\n")
    cat("=== Settings ===\n")
    cat("gene ID column:", if (!isTRUE(input$no_gene_id)) "yes (first column)" else "no (auto-assigned)", "\n")
    cat("analysis method:", input$analysis_method, "\n")
    cat("normalize method:", input$normalize_method, "\n")
    cat("log transform:", input$log_transform, "\n")
    cat("min count:", input$min_count, "\n")

    if (input$analysis_method == "auto") {
      cat("standard.min.group.n: ", input$standard_min_group_n, "\n")
    }

    if (input$analysis_method != "standard_anova") {
      cat("n.bin: ", input$n_bin, "\n")
      cat("spline df: ", input$df, "\n")
      cat("use_weighted_between: ", input$use_weighted_between, "\n")
      cat("trimming method: ", input$trim_method, "\n")
      if (input$trim_method == "dvalue") {
        cat("d.threshold: ", input$d_threshold, "\n")
      }
      cat("p-value method: ", input$p_method, "\n")
    }
  })
}

shinyApp(ui = ui, server = server)
