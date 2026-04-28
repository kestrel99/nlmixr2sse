#' Define an SSE alternative model
#'
#' `sseModel()` stores the model metadata that [runSSE()] will use when it
#' refits alternative hypotheses against each simulated dataset. It resolves
#' labels, estimation metadata, and fitted-versus-unfitted model inputs into a
#' common internal representation for the SSE workflow.
#'
#' @param model Required model specification. This may be either a fitted
#'   `nlmixr2FitCore` object or an unfitted model/UI object.
#' @param est Optional estimation method. Required when `model` is not a fitted
#'   object; otherwise defaults to the method recorded on the fit.
#' @param control Optional estimation control object. Required when `model` is
#'   not a fitted object; otherwise defaults to the control recorded on the fit.
#' @param label Optional short label used in SSE outputs. When `NULL`,
#'   [runSSE()] auto-generates `alt1`, `alt2`, and so on.
#' @param data Optional data override. Reserved for future SSE phases.
#'
#' @return An object of class `nlmixr2SSEModel`.
#' @export
sseModel <- function(
  model,
  est = NULL,
  control = NULL,
  label = NULL,
  data = NULL
) {
  if (missing(model)) {
    .abortSSE("{.arg model} is required.")
  }
  if (!is.null(label)) {
    .validateModelLabel(label)
  }

  is_fit <- inherits(model, "nlmixr2FitCore")
  if (is_fit) {
    .validateSSEFit(model, arg = "model")
    resolved_est <- if (is.null(est)) .fitField(model, "est") else est
    resolved_control <- if (is.null(control)) {
      .fitField(model, "control")
    } else {
      control
    }
    ui <- .getFitUi(model)
  } else {
    resolved_est <- est
    resolved_control <- control
    if (is.null(resolved_est)) {
      .abortSSE(
        "Supply {.arg est} when {.arg model} is not a fitted {.cls nlmixr2FitCore} object."
      )
    }
    if (is.null(resolved_control)) {
      .abortSSE(
        "Supply {.arg control} when {.arg model} is not a fitted {.cls nlmixr2FitCore} object."
      )
    }
    ui <- if (inherits(model, "rxUi")) {
      model
    } else {
      tryCatch(
        nlmixr2est::nlmixr2(model),
        error = function(e) {
          .abortSSE(
            "Could not parse {.arg model} into an rxUi object: {conditionMessage(e)}"
          )
        }
      )
    }
  }

  if (!.isScalarCharacter(as.character(resolved_est))) {
    .abortSSE("{.arg est} must resolve to a single estimation-method string.")
  }
  if (is.null(resolved_control)) {
    .abortSSE("{.arg control} must resolve to a non-NULL control object.")
  }

  structure(
    list(
      model = model,
      ui = ui,
      est = as.character(resolved_est),
      control = resolved_control,
      label = label,
      data = data,
      isFit = is_fit
    ),
    class = "nlmixr2SSEModel"
  )
}
