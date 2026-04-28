#' S3 methods for `nlmixr2SSE`
#'
#' `summary()` reshapes the tidy Phase 4 SSE outputs into wide,
#' human-readable tables. `print()` displays a compact run overview plus the
#' widened parameter and OFV summaries.
#'
#' @param object,x An `nlmixr2SSE` object.
#' @param ... Reserved for future method-specific arguments.
#'
#' @return `summary()` returns an object of class `summary.nlmixr2SSE`.
#'   `print()` methods return their input invisibly.
#' @name nlmixr2SSE-methods
NULL

#' @rdname nlmixr2SSE-methods
#' @export
summary.nlmixr2SSE <- function(object, ...) {
  structure(
    list(
      runInfo = object$runInfo,
      parameterSummary = .wideParameterSummary(object$parameterSummary),
      ofvSummary = .wideOfvSummary(object$ofvSummary),
      powerSummary = .wideOfvSummary(object$powerSummary),
      outputDir = object$outputDir
    ),
    class = "summary.nlmixr2SSE"
  )
}

#' @rdname nlmixr2SSE-methods
#' @export
print.summary.nlmixr2SSE <- function(x, ...) {
  cat("<summary.nlmixr2SSE>\n")
  if (is.list(x$runInfo) && length(x$runInfo) > 0L) {
    cat("  Model: ", x$runInfo$fitName, "\n", sep = "")
    cat("  Samples: ", x$runInfo$samples, "\n", sep = "")
    cat("  Parameter source: ", x$runInfo$parameterSource, "\n", sep = "")
    cat("  Output directory: ", x$outputDir, "\n", sep = "")
  }

  if (nrow(x$parameterSummary) > 0L) {
    cat("\nParameter summary\n")
    print(x$parameterSummary, row.names = FALSE)
  }

  if (nrow(x$ofvSummary) > 0L) {
    cat("\nOFV summary\n")
    print(x$ofvSummary, row.names = FALSE)
  }

  invisible(x)
}

#' @rdname nlmixr2SSE-methods
#' @export
print.nlmixr2SSE <- function(x, ...) {
  print(summary(x, ...))
  invisible(x)
}
