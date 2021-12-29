#' Calculate a 5 parameter logistic regression model with some start values from another model
#'
#' nls.multstart::nls_multstart is used to calculate a 5 parameter logistic regression.
#' Multiparametric regression requires to guess starting values. This guessing is done
#' based on standard values (aa, dd) and on a previous regression made by the bioplex
#' biomanager (cc, dd, gg).
#'
#' @param list one list entry as derived from plexr::read_bioplex
#' @param n_iter number of iterations for nls.multstart::nls_multstart
#' @param weights on which end of the regression does good fitting to standards
#' is required most (e.g. where do the measured values lie), lower range or upper range
#' or no weights (none); only used if FI is NULL.
#' @param FI provide custom weights instead of the easy-to-use weights parameter
#' @param model_name if NULL, just the model is returned; if model_name is not NULL
#' the model is added under that name to list and list is returned
#' @param ... additional parameters passed to nls.multstart::nls_multstart; e.g.
#' control = nls.lm.control(maxiter = 1024, maxfev=10000); or supp_errors = "Y"
#'
#' @return
#' @export
#'
#' @examples
alt_5pl_model <- function(list,
                          n_iter = 1000,
                          weights = c("lower", "upper", "none"),
                          FI = NULL,
                          model_name = NULL,
                          ...) {

  weights <- match.arg(weights, c("lower", "upper", "none"))

  start_lower <- c(aa = max(list$standard$FI),
                   dd = min(list$standard$FI)/1.5,
                   gg = list$std_curve_pars[["gg"]]/10,
                   bb = list$std_curve_pars[["bb"]]*10,
                   cc = list$std_curve_pars[["cc"]]/10)
  start_upper <- c(aa = 100*max(list$standard$FI),
                   dd = min(list$standard$FI)*1.5,
                   gg = list$std_curve_pars[["gg"]]*10,
                   bb = list$std_curve_pars[["bb"]]/10,
                   cc = list$std_curve_pars[["cc"]]*10)

  if (is.null(FI)) {
    if (weights == "lower") {
      FI <- rev(list$standard$FI)^2
    }
    if (weights == "upper") {
      FI <- list$standard$FI^2
    }
    if (weights == "none") {
      FI <- rep(1, length(list$standard$FI))
    }
  }


  model <- nls.multstart::nls_multstart(formula = FI ~ dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg,
                                        data = list$standard,
                                        iter = n_iter,
                                        start_lower = start_lower,
                                        start_upper = start_upper,
                                        modelweights = FI,
                                        ...)
  if (!is.null(model_name)) {
    list[[model_name]] <- model
    return(list)
  } else {
    return(model)
  }
}
