#' nlmixr2sse package
#'
#' Stochastic simulation and estimation tools for nlmixr2 fits.
#'
#' @keywords internal
"_PACKAGE"

# These are ggplot2 non-standard-evaluation column names referenced inside
# this package's plotting functions (plotSSEPpePower(), plotSSEPpeDiagnostics(),
# plotSSEParameterDraws(), etc.) -- mostly ggplot2::aes() mappings, but also a
# few (`delta`, `OFV`, `.`) from ggplot2::label_bquote() facet labels, which
# use the same NSE mechanism. None of these are undeclared global variables.
# R CMD check cannot tell NSE column references from real global-variable use
# and flags every one as "no visible binding for global variable"; this is the
# standard ggplot2 NSE workaround (see ?globalVariables). Every plotting
# function that adds a new aes()/label_bquote() column name should add it
# here too, alphabetically, or this check regresses.
utils::globalVariables(c(
  ".", "alpha", "comparison", "delta", "delta_ofv", "empirical", "estimate",
  "fitted", "lower", "model_label", "OFV", "power", "power_lower",
  "power_upper", "rate", "rate_lower", "rate_upper", "role", "study_size",
  "target_mean", "theoretical", "threshold", "truth", "upper", "value"
))
