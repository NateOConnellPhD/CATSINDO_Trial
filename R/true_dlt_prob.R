#' Generate a True DLT Probability Function
#'
#' Uses a logistic model fit to two known dose-toxicity points to define the true probability
#' of DLT as a function of threshold.
#'
#' @param t_low Numeric. Threshold where DLT probability is known (e.g., 0.10).
#' @param p_low Numeric. DLT probability at \code{t_low}.
#' @param t_target Numeric. Threshold assumed to be the MTD.
#' @param p_target Numeric. DLT probability at the MTD.
#'
#' @return A function \code{f(thresh)} returning the DLT probability.
#'
#' @export
true_dlt_prob <- function(thresh, dlt_model) {
  plogis(dlt_model$alpha + dlt_model$beta * thresh)
}
