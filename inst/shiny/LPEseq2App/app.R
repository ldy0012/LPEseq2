options(shiny.maxRequestSize = 100 * 1024^2)
# Show the actual error message (e.g. "Preprocessing failed: ...") instead of
# Shiny's generic sanitized message ("An error has occurred. Check your logs
# or contact the app author for clarification.").
options(shiny.sanitize.errors = FALSE)

library(shiny)
library(DT)
library(bslib)
library(LPEseq2)

# ------------------------------------------------------------
#    Theme: clean grayscale surface + single accent color (blue)
#    Fonts are loaded via Google's CDN at the CSS level (local = FALSE)
#    rather than downloaded/cached on the server, so app startup does
#    not depend on outbound network access from the R process.
# ------------------------------------------------------------
lpe_theme <- bslib::bs_theme(
  version = 5,
  bg = "#FAFAFA",
  fg = "#1F2328",
  primary = "#2563EB",
  secondary = "#6B7280",
  success = "#16A34A",
  danger = "#DC2626",
  warning = "#D97706",
  info = "#2563EB",
  base_font = bslib::font_google("Inter", local = FALSE),
  heading_font = bslib::font_google("Inter", wght = "600", local = FALSE),
  code_font = bslib::font_google("JetBrains Mono", local = FALSE),
  "border-radius" = "10px",
  "border-color" = "#E5E7EB"
)

ui <- fluidPage(
  theme = lpe_theme,
  tags$head(
    tags$style(HTML("
      body {
        color: #1F2328;
      }

      /* ---------- header ---------- */
      .lpe-header {
        padding: 20px 0 16px 0;
        margin-bottom: 24px;
        border-bottom: 1px solid #E5E7EB;
      }
      .lpe-header h1 {
        font-size: 1.6rem;
        font-weight: 700;
        margin: 0 0 4px 0;
        display: flex;
        align-items: center;
        gap: 10px;
      }
      .lpe-header h1 .fa, .lpe-header h1 svg {
        color: #2563EB;
      }
      .lpe-header p {
        margin: 0;
        color: #6B7280;
        font-size: 0.92rem;
      }

      /* ---------- sidebar card ---------- */
      .well {
        background-color: #FFFFFF;
        border: 1px solid #E5E7EB;
        border-radius: 12px;
        box-shadow: 0 1px 2px rgba(0,0,0,0.03);
        padding: 20px;
      }

      .lpe-section-title {
        display: flex;
        align-items: center;
        gap: 8px;
        font-size: 0.92rem;
        font-weight: 700;
        color: #1F2328;
        text-transform: uppercase;
        letter-spacing: 0.03em;
        margin: 24px 0 14px 0;
        padding-bottom: 8px;
        border-bottom: 2px solid rgba(37, 99, 235, 0.15);
      }
      .lpe-section-title:first-child { margin-top: 0; }
      .lpe-section-title .badge-num {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 20px;
        height: 20px;
        flex: 0 0 20px;
        border-radius: 50%;
        background: #2563EB;
        color: #FFFFFF;
        font-size: 0.7rem;
        font-weight: 700;
      }

      /* ---------- form controls ---------- */
      label, .control-label {
        font-weight: 600;
        font-size: 0.85rem;
        color: #374151;
      }

      .form-control, .selectize-input {
        border-radius: 8px;
        border-color: #D1D5DB;
      }
      .form-control:focus, .selectize-input.focus {
        border-color: #2563EB;
        box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.12);
      }

      .help-block {
        font-size: 0.78rem;
        color: #6B7280;
        line-height: 1.4;
      }

      /* ---------- buttons ---------- */
      #run {
        width: 100%;
        padding: 10px 0;
        font-weight: 700;
        font-size: 0.95rem;
        border-radius: 8px;
        border: none;
        background-color: #2563EB;
        box-shadow: 0 1px 2px rgba(37, 99, 235, 0.3);
      }
      #run:hover, #run:focus {
        background-color: #1D4ED8;
      }

      #download_results {
        width: 100%;
        border-radius: 8px;
        font-weight: 600;
        background-color: #FFFFFF;
        border: 1px solid #D1D5DB;
        color: #1F2328;
      }
      #download_results:hover {
        background-color: #F3F4F6;
        border-color: #9CA3AF;
      }

      /* ---------- run status ---------- */
      .run-status-running {
        color: #2563EB;
        font-weight: 700;
        padding: 8px 0;
      }
      .run-status-done {
        color: #16A34A;
        font-weight: 700;
        padding: 8px 0;
      }
      .run-status-error {
        color: #DC2626;
        font-weight: 700;
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

      /* ---------- tabs / main panel ---------- */
      .nav-tabs {
        border-bottom: 1px solid #E5E7EB;
        gap: 2px;
        flex-wrap: wrap;
      }
      .nav-tabs .nav-link {
        border: none;
        border-radius: 8px 8px 0 0;
        color: #6B7280;
        font-weight: 600;
        font-size: 0.84rem;
        padding: 10px 14px;
      }
      .nav-tabs .nav-link.active {
        color: #2563EB;
        background-color: #FFFFFF;
        border-bottom: 2px solid #2563EB;
      }
      .nav-tabs .nav-link:hover:not(.active) {
        color: #1F2328;
        background-color: #F3F4F6;
        border-color: transparent;
      }

      .tab-content {
        background-color: #FFFFFF;
        border: 1px solid #E5E7EB;
        border-top: none;
        border-radius: 0 0 12px 12px;
        padding: 22px;
      }

      /* ---------- misc ---------- */
      hr {
        border-top: 1px solid #E5E7EB;
        margin: 20px 0;
      }

      .gene-id-note {
        background-color: #FEF9C3;
        border-left: 4px solid #EAB308;
        border-radius: 8px;
        padding: 10px 14px;
        margin-top: 6px;
        font-size: 0.85rem;
      }

      pre, .shiny-text-output {
        border-radius: 8px;
        background-color: #F9FAFB;
        border: 1px solid #E5E7EB;
      }

      table.dataTable {
        font-size: 0.85rem;
      }

      /* ---------- top control bar ---------- */
      .lpe-topbar {
        margin-bottom: 20px;
      }
      .lpe-topbar .well {
        margin-bottom: 12px;
      }
      .accordion-button:not(.collapsed) {
        color: #2563EB;
        background-color: rgba(37, 99, 235, 0.06);
      }
      .accordion-button:focus {
        box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.12);
      }
    "))
  ),

  div(
    class = "lpe-header",
    h1(icon("dna"), "LPEseq2"),
    p("Local Pooled Error-Based ANOVA for RNA-Seq Count Data")
  ),

  div(
    class = "lpe-topbar",
    wellPanel(
      fluidRow(
        column(
          3,
          fileInput("counts_file", "Upload counts file", accept = c(".csv", ".tsv", ".txt")),
          checkboxInput(
            "no_gene_id",
            "First column is NOT a gene identifier (auto-assign gene IDs)",
            value = FALSE
          ),
          uiOutput("gene_id_warning_ui")
        ),

        column(
          3,
          fileInput("meta_file", "Upload metadata file", accept = c(".csv", ".tsv", ".txt")),
          uiOutput("group_var_ui")
        ),

        column(
          2,
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
          selectInput(
            "normalize_method",
            "Normalization method",
            choices = c("library_size", "TMM", "DESeq2", "none"),
            selected = "TMM"
          )
        ),

        column(
          2,
          selectInput(
            "variance_eval",
            "Variance evaluation method",
            choices = c("Grand mean (default)" = "grand_mean",
                        "Per-group (Welch)" = "per_group"),
            selected = "grand_mean"
          ),
          checkboxInput(
            "log_transform",
            "Log2 transform",
            value = TRUE
          )
        ),

        column(
          2,
          actionButton(
            "run",
            "Run Analysis",
            icon = icon("play"),
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
        )
      )
    ),

    bslib::accordion(
      id = "advanced_options_accordion",
      open = FALSE,
      bslib::accordion_panel(
        title = "Advanced options",
        icon = icon("sliders-h"),

        fluidRow(
          column(
            3,
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
            )
          ),

          column(
            3,
            helpText(
              "Per-group evaluates variance separately for each group's own mean ",
              "expression level, matching LPEseq1's approach for 2-group comparisons. ",
              "Recommended when group means differ substantially in intensity."
            )
          ),

          column(
            2,
            numericInput(
              "min_count",
              "Minimum count",
              value = 5,
              min = 0
            )
          ),

          column(
            2,
            numericInput(
              "prior_count",
              "Pseudo count",
              value = 1,
              min = 0
            )
          )
        ),

        conditionalPanel(
          condition = "input.analysis_method != 'standard_anova'",

          tags$hr(),

          fluidRow(
            column(
              2,
              numericInput(
                "n_bin",
                "Number of bins",
                value = 100,
                min = 5
              )
            ),

            column(
              2,
              numericInput(
                "df",
                "Spline degrees of freedom",
                value = 10,
                min = 2
              )
            ),

            column(
              3,
              selectInput(
                "trim_method",
                "Pairwise outlier trimming method",
                choices = c(
                  "Pooled bin-wise IQR" = "iqr",
                  "Fixed D-value threshold (LPEseq1)" = "dvalue",
                  "None" = "none"
                ),
                selected = "dvalue"
              ),

              helpText(
                "The IQR method pools within-group and between-group-derived pairwise values, ",
                "divides them into expression-intensity A-bins, and applies the conventional 1.5 within each bin. ",
                "Outlier detection is performed on the M-value scale used for variance trend estimation."
              )
            ),

            column(
              2,
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
                  "Applies a fixed threshold on the M scale (the rescaled value actually ",
                  "used for variance estimation), converted from this D-value setting as ",
                  "threshold/sqrt(2), as in LPEseq1's non-replicate outlier procedure ",
                  "(Gim et al. 2016). Any pairwise value whose |M| exceeds the converted ",
                  "threshold is removed, regardless of expression-intensity bin or group ",
                  "size. Default 1.2 was empirically tuned on specific benchmark datasets; ",
                  "consider adjusting for your data."
                )
              )
            ),

            column(
              3,
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
            )
          )
        )
      )
    )
  ),

  div(
    style = "margin-top: 20px;",
    tabsetPanel(
      tabPanel(
        "Instructions",
        icon = icon("info-circle"),

        h4(icon("rocket"), " Quick start"),
        tags$ol(
          tags$li("Upload your ", strong("counts file"), " (genes as rows, samples as columns)."),
          tags$li("Upload your ", strong("metadata file"), " (samples as rows, variables such as group as columns)."),
          tags$li("Select the ", strong("group variable"), " you want to compare (e.g. \"group\")."),
          tags$li("Adjust the analysis options if needed — the defaults work for most datasets."),
          tags$li("Click ", strong("Run Analysis"), " and check the ", strong("Results"), " tab once it finishes.")
        ),

        tags$hr(),

        h4(icon("file-alt"), " Input file requirements"),
        tags$ul(
          tags$li(strong("Counts file: "), "genes as rows, samples as columns."),
          tags$li(strong("Metadata file: "), "samples as rows, variables (e.g. group) as columns."),
          tags$li("The ", strong("column names"), " of the counts file must exactly match the ", strong("row names"), " of the metadata file (i.e. the sample names)."),
          tags$li("If your counts file has no gene identifier column, check ", em("“First column is NOT a gene identifier”"), " in the sidebar to auto-assign gene IDs (gene_1, gene_2, ...).")
        ),

        tags$hr(),

        h4(icon("table"), " Counts file format"),
        tags$ul(
          tags$li("Accepted formats: ", code(".csv"), ", ", code(".tsv"), ", ", code(".txt"),
                  " (comma, tab, semicolon, or pipe separated — detected automatically)."),
          tags$li("First column: gene identifiers (used as row names)."),
          tags$li("Remaining columns: one column per sample, integer counts recommended."),
          tags$li("First row: header with sample names. No missing values allowed.")
        ),
        p(strong("Example:")),
        div(
          style = "overflow-x: auto;",
          tags$table(
            class = "table table-sm table-bordered",
            style = "max-width: 480px;",
            tags$thead(
              tags$tr(
                tags$th("gene"), tags$th("sample1"), tags$th("sample2"),
                tags$th("sample3"), tags$th("sample4")
              )
            ),
            tags$tbody(
              tags$tr(tags$td("gene1"), tags$td("100"), tags$td("120"), tags$td("80"), tags$td("95")),
              tags$tr(tags$td("gene2"), tags$td("50"), tags$td("60"), tags$td("55"), tags$td("70")),
              tags$tr(tags$td("gene3"), tags$td("10"), tags$td("15"), tags$td("30"), tags$td("28"))
            )
          )
        ),

        h4(icon("list"), " Metadata file format"),
        tags$ul(
          tags$li("Accepted formats: ", code(".csv"), ", ", code(".tsv"), ", ", code(".txt"),
                  " (comma, tab, semicolon, or pipe separated — detected automatically)."),
          tags$li("First column: sample names (must match the column names of the counts file)."),
          tags$li("Remaining columns: one column per variable (e.g. group, batch)."),
          tags$li("First row: header with variable names.")
        ),
        p(strong("Example:")),
        div(
          style = "overflow-x: auto;",
          tags$table(
            class = "table table-sm table-bordered",
            style = "max-width: 320px;",
            tags$thead(
              tags$tr(tags$th("sample"), tags$th("group"))
            ),
            tags$tbody(
              tags$tr(tags$td("sample1"), tags$td("Control")),
              tags$tr(tags$td("sample2"), tags$td("Control")),
              tags$tr(tags$td("sample3"), tags$td("Treatment")),
              tags$tr(tags$td("sample4"), tags$td("Treatment"))
            )
          )
        )
      ),

      tabPanel(
        "Counts preview",
        icon = icon("table"),
        DTOutput("counts_preview")
      ),

      tabPanel(
        "Metadata preview",
        icon = icon("list"),
        DTOutput("meta_preview")
      ),

      tabPanel(
        "Results",
        icon = icon("chart-bar"),
        DTOutput("results_table")
      ),

      tabPanel(
        "Method info",
        icon = icon("cog"),
        verbatimTextOutput("method_info")
      ),

      tabPanel(
        "Variance trend info",
        icon = icon("chart-line"),
        verbatimTextOutput("trend_info"),
        DTOutput("base_var_table")
      ),

      tabPanel(
        "Trimming info",
        icon = icon("filter"),
        verbatimTextOutput("trim_info"),
        DTOutput("trim_table")
      ),

      tabPanel(
        "Spline Plot",
        icon = icon("chart-area"),
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
        icon = icon("terminal"),
        verbatimTextOutput("log_text"),
      )
    )
  )
)

server <- function(input, output, session) {

  run_state <- reactiveVal("idle")  # idle / running / done / error
  gene_id_warning <- reactiveVal(NULL)

  # Holds the outcome of the most recent Run click: a list with $status
  # ("done" / "error" / "idle") and either $data (the result data.frame, on
  # success) or $condition (the original R/validation condition, on
  # failure). The analysis is computed inside a plain observeEvent() rather
  # than a lazy eventReactive(): observers are never suspended just because
  # the tab that would display their result isn't currently visible, so
  # this guarantees the analysis always runs the instant "Run Analysis" is
  # clicked, no matter which tab the user is on at the time.
  analysis_state <- reactiveVal(list(status = "idle"))

  observeEvent(input$run, {
    run_state("running")

    outcome <- tryCatch({

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
      lpe_trim_method <- if (is.null(input$trim_method)) "dvalue" else input$trim_method
      lpe_use_weighted_between <- if (is.null(input$use_weighted_between)) TRUE else input$use_weighted_between
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
        p.method           = lpe_p_method,
        variance.eval = input$variance_eval
      )

      list(status = "done", data = res)
    }, error = function(e) {
      if (inherits(e, "shiny.silent.error")) {
        list(status = "idle", condition = e)
      } else {
        list(status = "error", condition = e)
      }
    })

    analysis_state(outcome)
    run_state(outcome$status)
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
      class = "gene-id-note",
      icon("exclamation-triangle"),
      strong(" Note: "),
      msg
    )
  })

  analysis_result <- reactive({
    state <- analysis_state()
    if (!is.null(state$condition)) {
      stop(state$condition)
    }
    req(identical(state$status, "done"))
    state$data
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
  #   xlab_txt <- if (attr(res, "variance.eval") == "per_group") "log2(Welch T-statistic)" else "log2(MS_between)"
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
    cat("   Note: Counts columns must match metadata row names.\n")
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

  # ------------------------------------------------------------
  # By default, Shiny suspends (pauses/cancels) computation for
  # outputs on tabs that aren't currently visible. If the user
  # switches to a results tab while LPE_ANOVA() is still running
  # in the background, that suspend/resume cycle can interrupt
  # the in-progress computation and surface an error. Disabling
  # suspendWhenHidden for every output that depends on
  # analysis_result() keeps them "live" regardless of which tab
  # is showing, so switching tabs mid-run no longer interrupts
  # the analysis.
  # ------------------------------------------------------------
  outputOptions(output, "results_table", suspendWhenHidden = FALSE)
  outputOptions(output, "method_info", suspendWhenHidden = FALSE)
  outputOptions(output, "trend_info", suspendWhenHidden = FALSE)
  outputOptions(output, "base_var_table", suspendWhenHidden = FALSE)
  outputOptions(output, "trim_info", suspendWhenHidden = FALSE)
  outputOptions(output, "trim_table", suspendWhenHidden = FALSE)
  outputOptions(output, "spline_plot", suspendWhenHidden = FALSE)
}

shinyApp(ui = ui, server = server)
