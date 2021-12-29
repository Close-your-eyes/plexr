#' 5-parameter logistic regression prediction
#'
#' Alternative to stats::predict. Provide concentration (Conc) or Fluorescence intensity (FI)
#' to have the other one returned according to the model.
#'
#' If FI is provided, dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg is used;
#' if Conc is provided dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg is used.
#' Alternatively stats::predict(nls_model, newdata = data.frame(Conc = 1:10))
#'
#' @param nls_model non-linear least squares regression model
#' @param Conc concentration values
#' @param FI fluorescence intensities
#'
#' @return
#' @export
#'
#' @examples
five_par_log_regress <- function(nls_model, Conc = NULL, FI = NULL) {

  if (missing(nls_model)) {
    stop("Please provide a model.")
  }

  if (!is.null(FI) && !is.null(Conc)) {
    stop("Please provide Conc or FI, but not both.")
  }

  model <- as.data.frame(broom::tidy(nls_model))
  dd <- model[which(model$term == "dd"),"estimate"]
  aa <- model[which(model$term == "aa"),"estimate"]
  cc <- model[which(model$term == "cc"),"estimate"]
  bb <- model[which(model$term == "bb"),"estimate"]
  gg <- model[which(model$term == "gg"),"estimate"]

  if (any(lengths(list(dd, aa, cc, bb, gg)) == 0)) {
    stop("At least one parameter is missing in the model.")
  }

  if (is.null(FI)) {
    return(dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg)
  }
  if (is.null(Conc)) {
    return( ((((aa - dd) / (FI + dd))^(1/gg) - 1)^(1/bb))*cc )
  }
}

