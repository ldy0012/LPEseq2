fixbounds.predict.smooth.spline <- function(object, x) {

  if (!inherits(object, "smooth.spline")) {
    stop("object must be a smooth.spline object")
  }

  if (is.null(object$x)) {
    stop("object must be a smooth.spline object with x values")
  }

  if (!is.numeric(x)) {
    stop("x must be numeric")
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop("x contains NA, NaN, or infinite values")
  }

  valid_x <- object$x[is.finite(object$x)]

  if (length(valid_x) == 0) {
    stop("object$x contains no finite values")
  }

  xmin <- min(valid_x)
  xmax <- max(valid_x)

  x2 <- pmin(pmax(x, xmin), xmax)

  stats::predict(object, x2)
}
