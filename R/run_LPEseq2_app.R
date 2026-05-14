#' Launch the LPEseq2 Shiny App
#'
#' This function launches the interactive Shiny web application for LPEseq2.
#'
#' @return No return value. Opens a Shiny application.
#'
#' @export
run_LPEseq2_app <- function() {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required to run the LPEseq2 app.")
  }
  
  app_dir <- system.file("shiny", "LPEseq2App", package = "LPEseq2")
  
  if (app_dir == "") {
    stop("Could not find the LPEseq2 Shiny app directory.")
  }
  
  shiny::runApp(app_dir)
}