#' Create a True Dose-Toxicity Model Using a Logistic or Probit Link
#'
#' Constructs a dose-toxicity function passing through two known dose-DLT pairs — the anchor point and the assumed true best threshold at the target DLT rate.
#' Output list contains all information about the true dose-toxicity response function assumed for simulations.
#' Supports either a logit or probit link function.
#'
#' @param anchor_thresh Numeric. Anchor dose (e.g., 0.10).
#' @param anchor_dlt_rate Numeric. DLT probability at the anchor dose (e.g., 0.10).
#' @param true_mtd Numeric. Threshold corresponding to the target DLT probability (e.g., 0.18).
#' @param true_dlt_rate Numeric. Target DLT probability at the true MTD (e.g., 0.15).
#' @param max_thresh Numeric. Maximum threshold value to plot.
#' @param link Character. Link function to use: "logit" or "probit".
#' @param plot Logical. Whether to plot the resulting curve.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{true_prob_dlt}}{Function mapping thresholds to DLT probabilities.}
#'   \item{\code{alpha}}{Intercept of the model.}
#'   \item{\code{beta}}{Slope of the model.}
#'   \item{\code{anchor_thresh}}{The anchor threshold.}
#'   \item{\code{anchor_dlt_rate}}{DLT probability at the anchor threshold.}
#'   \item{\code{true_mtd}}{The dose at the target DLT probability (true MTD).}
#'   \item{\code{true_dlt_rate}}{Target DLT probability.}
#'   \item{\code{target_dlt_rate}}{Alias for \code{true_dlt_rate}.}
#'   \item{\code{max_thresh}}{Maximum threshold used in plotting.}
#'   \item{\code{link}}{Link function used: "logit" or "probit".}
#' }
#'
#' @export
true_model_old <- function(
    anchor_thresh,
    anchor_dlt_rate,
    true_mtd,
    true_dlt_rate,
    max_thresh,
    link = c("logit", "probit"),
    plot = TRUE
) {
  link <- match.arg(link)
  
  # Choose link functions
  qlink <- if (link == "logit") qlogis else qnorm
  plink <- if (link == "logit") plogis else pnorm
  
  # Convert probabilities to link scale
  link_anchor <- qlink(anchor_dlt_rate)
  link_target <- qlink(true_dlt_rate)
  
  # Estimate slope and intercept
  beta <- (link_target - link_anchor) / (true_mtd - anchor_thresh)
  alpha <- link_anchor - beta * anchor_thresh
  
  # Define dose-toxicity function
  true_prob_dlt <- function(thresh) plink(alpha + beta * thresh)
  
  # Optional plot
  if (plot) {
    threshold_range <- seq(min(anchor_thresh, true_mtd) - 0.02, max_thresh, by = 0.001)
    
    plot(threshold_range, true_prob_dlt(threshold_range),
         type = "l", lwd = 2, col = "blue",
         ylab = "DLT Probability", xlab = "Threshold",
         main = paste("True Dose-Toxicity Curve (", link, " link)", sep = ""))
    
    abline(h = true_dlt_rate, v = true_mtd, col = "red", lty = 2)
    
    text(true_mtd, true_dlt_rate - 0.01,
         labels = paste0("True DLT Rate = ", true_dlt_rate, "\nTrue MTD = ", true_mtd),
         col = "red", pos = 4)
    
    legend("topleft",
           legend = c("True DLT Curve", "Target DLT / True MTD"),
           col = c("blue", "red"),
           lty = c(1, 2), lwd = c(2, 1))
  }
  
  list(
    true_prob_dlt = true_prob_dlt,
    alpha = alpha,
    beta = beta,
    anchor_thresh = anchor_thresh,
    anchor_dlt_rate = anchor_dlt_rate,
    true_mtd = true_mtd,
    true_dlt_rate = true_dlt_rate,
    max_thresh = max_thresh,
    link = link
  )
}

