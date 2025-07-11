#' Construct prior information for dose-escalation model
#'
#' Generates prior parameters for a Bayesian dose–DLT logistic or probit regression model.
#' Allows dose rescaling for numerical stability and optionally forces the prior curve 
#' through a known MTD and DLT rate.
#'
#' @param doses Numeric vector of raw dose levels used to define the prior.
#' @param probs Numeric vector of prior DLT probabilities corresponding to each dose.
#' @param true_mtd The known or hypothesized MTD dose (on the raw scale).
#' @param true_dlt_rate The DLT probability at the true MTD (e.g., 0.15).
#' @param force_through_mtd Logical. If TRUE, forces the prior curve through the true MTD and DLT rate.
#' @param int_prior Standard deviation of the prior on the intercept (alpha).
#' @param slope_prior Standard deviation of the prior on the slope (beta). If NULL, defaults to 1 / range of scaled doses.
#' @param link Link function to use. Either "logit" or "probit".
#' @param dose_scale_factor Multiplicative factor applied after subtracting the minimum dose, to stabilize estimation (e.g., 100).
#' @param plot Logical. If TRUE, plots the prior curve with uncertainty bands.
#' @param plot_range Range of doses over which to plot the prior (default: range of input doses).
#' @param n_draws Number of samples drawn from the prior to compute credible bands in the plot.
#'
#' @return A list with prior parameters and helper functions:
#' \describe{
#'   \item{true_prob_dlt}{A function that maps raw doses to DLT probabilities under the prior}
#'   \item{alpha}{Intercept of the scaled logistic/probit model}
#'   \item{beta}{Slope of the scaled model}
#'   \item{int_prior}{Intercept prior SD}
#'   \item{slope_prior}{Slope prior SD}
#'   \item{dose_base}{Minimum dose used to scale}
#'   \item{doses}{Input raw dose vector}
#'   \item{probs}{Input DLT probabilities}
#'   \item{true_mtd}{Supplied true MTD}
#'   \item{true_dlt_rate}{DLT rate at the true MTD}
#'   \item{dose_scale_factor}{Scaling multiplier}
#'   \item{link}{Link function used ("logit" or "probit")}
#' }
#'
#' @examples
#' prior_info <- get_prior(
#'   doses = c(0.1, 0.13, 0.15, 0.17, 0.2),
#'   probs = c(0.1, 0.12, 0.15, 0.23, 0.3),
#'   true_mtd = 0.15,
#'   true_dlt_rate = 0.15,
#'   force_through_mtd = TRUE,
#'   int_prior = 0.3,
#'   slope_prior = 1,
#'   dose_scale_factor = 100,
#'   link = "logit",
#'   plot = TRUE
#' )
#'
#' # Evaluate prior at specific dose
#' prior_info$true_prob_dlt(0.15)
#'
#' @export

get_prior <- function(
    doses,
    probs,
    true_mtd,
    true_dlt_rate,
    force_through_mtd = FALSE,
    int_prior         = 0.3,
    slope_prior       = NULL,
    link              = c("logit", "probit"),
    dose_scale_factor = 100,
    plot              = TRUE,
    plot_range        = range(doses),
    n_draws           = 1000
) {
  stopifnot(length(doses) == length(probs))
  stopifnot(all(probs > 0 & probs < 1))
  
  link   <- match.arg(link)
  qlink  <- if (link == "logit") qlogis else qnorm
  plink  <- if (link == "logit") plogis else pnorm
  
  # Scale doses: shift + rescale
  dose_base <- min(doses)
  sdoses    <- (doses    - dose_base) * dose_scale_factor
  smtd      <- (true_mtd - dose_base) * dose_scale_factor
  
  # Auto-scale slope prior if not provided
  range_sdose <- diff(range(sdoses))
  if (is.null(slope_prior)) {
    slope_prior <- 1 / range_sdose
  }
  
  if (force_through_mtd) {
    target_link <- qlink(true_dlt_rate)
    loss_fn <- function(beta) {
      alpha <- target_link - beta * smtd
      pred  <- plink(alpha + beta * sdoses)
      sum((pred - probs)^2)
    }
    opt   <- optimize(loss_fn, interval = c(-100, 100))
    beta  <- opt$minimum
    alpha <- target_link - beta * smtd
  } else {
    logits <- qlink(probs)
    fit    <- lm(logits ~ sdoses)
    alpha  <- coef(fit)[1]
    beta   <- coef(fit)[2]
  }
  
  # Define prior predictive function
  true_prob_dlt <- function(dose_raw) {
    dose_scaled <- (dose_raw - dose_base) * dose_scale_factor
    plink(alpha + beta * dose_scaled)
  }
  
  # Plot
  if (plot) {
    rr  <- seq(plot_range[1], plot_range[2], length.out = 1000)
    srr <- (rr - dose_base) * dose_scale_factor
    
    i_samp <- rnorm(n_draws, mean = alpha, sd = int_prior)
    b_samp <- rnorm(n_draws, mean = beta,  sd = slope_prior)
    
    m <- sapply(seq_len(n_draws), function(i) plink(i_samp[i] + b_samp[i] * srr))
    ci_df <- tibble::tibble(
      dose  = rr,
      mean  = rowMeans(m),
      lower = apply(m, 1, quantile, 0.025),
      upper = apply(m, 1, quantile, 0.975)
    )
    
    plot(ci_df$dose, ci_df$mean, type = "l", lwd = 2, col = "darkred",
         xlab = "Dose", ylab = "DLT Probability",
         main = paste0("Prior Curve (", link, ")"),
         ylim = c(0, min(1, max(ci_df$upper) + 0.05)))
    polygon(c(ci_df$dose, rev(ci_df$dose)),
            c(ci_df$lower, rev(ci_df$upper)),
            col = adjustcolor("darkred", alpha.f = 0.15), border = NA)
    points(doses, probs, pch = 19, col = "red")
    abline(h = true_dlt_rate, v = true_mtd, col = "gray40", lty = 2)
  }
  
  list(
    true_prob_dlt     = true_prob_dlt,
    alpha             = alpha,
    beta              = beta,
    int_prior         = int_prior,
    slope_prior       = slope_prior,
    dose_base         = dose_base,
    doses             = doses,
    probs             = probs,
    true_mtd          = true_mtd,
    true_dlt_rate     = true_dlt_rate,
    dose_scale_factor = dose_scale_factor,
    link              = link
  )
}
