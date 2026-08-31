#' nlmixr2sse package
#'
#' Stochastic simulation and estimation tools for nlmixr2 fits.
#'
#' @keywords internal
"_PACKAGE"

# These are ggplot2::aes() column names referenced by non-standard
# evaluation inside this package's plotting functions (plotSSEPpePower(),
# plotSSEPpeDiagnostics(), plotSSEParameterDraws(), etc.) -- not undeclared
# global variables. R CMD check cannot tell the difference and flags every
# one as "no visible binding for global variable"; this is the standard
# ggplot2 NSE workaround (see ?globalVariables). Every plotting function that
# adds a new aes() column name should add it here too, alphabetically, or
# this check regresses.
utils::globalVariables(c(
  ".", "alpha", "comparison", "delta", "delta_ofv", "empirical", "estimate",
  "fitted", "lower", "model_label", "OFV", "power", "power_lower",
  "power_upper", "rate", "rate_lower", "rate_upper", "role", "study_size",
  "target_mean", "theoretical", "threshold", "truth", "upper", "value"
))
