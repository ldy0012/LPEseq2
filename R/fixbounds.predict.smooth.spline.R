fixbounds.predict.smooth.spline <- function(object, x) {

  if (is.null(object$x)) {
    stop("object must be a smooth.spline object with x values")
  }

  if (!is.numeric(x)) {
    stop("x must be numeric")
  }

  xmin <- min(object$x, na.rm = TRUE)
  xmax <- max(object$x, na.rm = TRUE)

  x2 <- pmin(pmax(x, xmin), xmax)

  stats::predict(object, x2)
}
