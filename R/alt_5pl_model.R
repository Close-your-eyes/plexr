#' Calculate a guided 5 parameter logistic regression model with some start values from another model
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
#' @param ... additional parameters passed to nls.multstart::nls_multstart; e.g. or supp_errors = "Y"
#' @param upper_limit_dd limit dd to the x-th smallest value of y (FI); dd is the lower limit of
#' values that can be  supplied to the inverse function of the model; so a low dd allows to supply lower
#' x (Conc) values to convert to FI
#' @param exag_factor factor used to increase range of allowed range of gg, bb, cc from std_curve_pars
#'
#' @return plain nls_model when is.null(model_name); or appended list with nls_model named model_name
#' @export
#'
#' @examples
alt_5pl_model <- function(list,
                          n_iter = 1000,
                          weights = c("lower", "upper", "none"),
                          model_name = NULL,
                          exag_factor = 10,
                          upper_limit_dd = NULL,
                          ...) {

  if (!is.numeric(weights)) {
    weights <- match.arg(weights, c("lower", "upper", "none"))
  }

  # just printing to find problematic analytes easier from outside loop
  if ("standard" %in% names(list)) {
    if ("sheet" %in% names(list[["standard"]])) {
      print(unique(list[["standard"]][["sheet"]]))
    } else if ("analyte" %in% names(list[["standard"]])) {
      print(unique(list[["standard"]][["analyte"]]))
    }
  } else if("samples" %in% names(list)) {
    if ("sheet" %in% names(list[["samples"]])) {
      print(unique(list[["samples"]][["sheet"]]))
    } else if ("analyte" %in% names(list[["samples"]])) {
      print(unique(list[["samples"]][["analyte"]]))
    }
  }

  start_lower <- c(aa = max(list$standard$FI),
                   dd = min(list$standard$FI)/1.5,
                   #gg = 0.1,
                   #bb = -2,
                   #cc = 100
                   gg = list$std_curve_pars[["gg"]]/exag_factor,
                   bb = list$std_curve_pars[["bb"]]*exag_factor,
                   cc = list$std_curve_pars[["cc"]]/exag_factor
  )
  start_upper <- c(aa = 100*max(list$standard$FI),
                   dd = min(list$standard$FI)*1.5,
                   #gg = 5,
                   #bb = -0.1,
                   #cc = 10000
                   gg = list$std_curve_pars[["gg"]]*exag_factor,
                   bb = list$std_curve_pars[["bb"]]/exag_factor,
                   cc = list$std_curve_pars[["cc"]]*exag_factor
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

  upper <- stats::setNames(c(rep(Inf, length(start_upper))), nm = names(start_upper))
  if (!is.null(upper_limit_dd)) {
    if (!is.numeric(upper_limit_dd) || length(upper_limit_dd) > 1) {
      stop("upper_limit_dd has to be a numeric (integer) of length 1.")
    }
    upper[["dd"]] <- sort(list$standard$FI)[upper_limit_dd]
  }

  nls_model <- nls.multstart::nls_multstart(formula = FI ~ dd + (aa - dd) / ((1 + (Conc / cc)^bb))^gg,
                                            data = list$standard,
                                            iter = n_iter,
                                            start_lower = start_lower,
                                            start_upper = start_upper,
                                            upper = upper,
                                            modelweights = w,
                                            control = minpack.lm::nls.lm.control(maxiter = 1024, maxfev = 10000), ...)

  model <- as.data.frame(broom::tidy(nls_model))
  dd <- model[which(model$term == "dd"),"estimate"]
  aa <- model[which(model$term == "aa"),"estimate"]
  cc <- model[which(model$term == "cc"),"estimate"]
  bb <- model[which(model$term == "bb"),"estimate"]
  gg <- model[which(model$term == "gg"),"estimate"]

  if (dd > aa) {
    warning("dd is greater than aa which is unexpected/should not be and may cause problems later on.")
  }
  if (cc < start_lower[["cc"]] || cc > start_upper[["cc"]]) {
    print("cc outside of start_lower <-> start_upper range.")
  }
  if (bb < start_lower[["bb"]] || bb > start_upper[["bb"]]) {
    print("bb outside of start_lower <-> start_upper range.")
  }
  if (gg < start_lower[["gg"]] || gg > start_upper[["gg"]]) {
    print("gg outside of start_lower <-> start_upper range.")
  }

  if (!is.null(model_name)) {
    list[[model_name]] <- nls_model
    return(list)
  } else {
    return(nls_model)
  }
}
