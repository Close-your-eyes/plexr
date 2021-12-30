#' Calculate a 5 parameter logistic regression model with some start values from another model
#'
#' nls.multstart::nls_multstart is used to calculate a 5 parameter logistic regression.
#' Multiparametric regression requires to guess starting values. This guessing is done
#' based on standard values (for aa, dd) and on a previous regression made by the bioplex
#' biomanager (for cc, dd, gg). Sometimes, results of model fitting have not been
#' consistent. Try to vary n_iter or simply repeat the fitting.
#'
#' Model for fitting is FI ~ dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg.
#' dd corresponds to the horizontal asymptote for x -> -Inf,
#' aa corresponds to the horizontal asymptote for x -> Inf.
#' (https://laustep.github.io/stlahblog/posts/5pl.html)
#'
#' @param list one list entry as derived from plexr::read_bioplex
#' @param n_iter number of iterations for nls.multstart::nls_multstart
#' @param weights on which end of the regression does good fitting to standards
#' is required most (e.g. where do the measured values lie), lower range or upper range
#' or no weights (none); alternatively provide a numeric vector of weights -
#' has to be same length as list$standard$FI; rows of list$standard will be order increasingly
#' by FI before fitting; take this into account when providing numeric weights.
#' @param model_name if NULL, just the model is returned; if model_name is not NULL
#' the model is added under that name to list and list is returned
#' @param ... additional parameters passed to nls.multstart::nls_multstart; e.g.
#' control = minpack.lm::nls.lm.control(maxiter = 1024, maxfev=10000); or supp_errors = "Y"
#'
#' @return plain nls_model when is.null(model_name); or appended list with nls_model named model_name
#' @export
#'
#' @examples
alt_5pl_model <- function(list,
                          n_iter = 1000,
                          weights = c("lower", "upper", "none"),
                          model_name = NULL,
                          ...) {

  if (!is.numeric(weights)) {
    weights <- match.arg(weights, c("lower", "upper", "none"))
  }

  start_lower <- c(aa = max(list$standard$FI),
                   dd = min(list$standard$FI)/1.5,
                   #gg = 0.1,
                   #bb = -2,
                   #cc = 100
                   gg = list$std_curve_pars[["gg"]]/10,
                   bb = list$std_curve_pars[["bb"]]*10,
                   cc = list$std_curve_pars[["cc"]]/10
  )
  start_upper <- c(aa = 100*max(list$standard$FI),
                   dd = min(list$standard$FI)*1.5,
                   #gg = 5,
                   #bb = -0.1,
                   #cc = 10000
                   gg = list$std_curve_pars[["gg"]]*10,
                   bb = list$std_curve_pars[["bb"]]/10,
                   cc = list$std_curve_pars[["cc"]]*10
  )

  list$standard <- list$standard[order(list$standard$FI),]

  if (!is.numeric(weights)) {
    if (weights == "lower") {
      list$standard$w <- rev(list$standard$FI)^2
    } else if (weights == "upper") {
      list$standard$w <- list$standard$FI^2
    } else if (weights == "none") {
      list$standard$w <- rep(1, times = length(list$standard$FI))
    }
  } else {
    if (length(weights) != length(list$standard$FI)) {
      stop("weights has to have the same length as list$standard$FI.")
    }
    list$standard$w <- weights
  }

  model <- nls.multstart::nls_multstart(formula = FI ~ dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg,
                                        data = list$standard,
                                        iter = n_iter,
                                        start_lower = start_lower,
                                        start_upper = start_upper,
                                        modelweights = w,
                                        ...)


  if (!is.null(model_name)) {
    list[[model_name]] <- model
    return(list)
  } else {
    return(model)
  }
}
