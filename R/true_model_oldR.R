
#' Fit a Logistic or Probit Dose-Toxicity Curve Constrained to Pass Through a Given MTD Point
#'
#' Fits a dose-toxicity curve using a logistic or probit link, constrained to pass
#' through a user-specified MTD (dose, DLT rate) pair. Minimizes squared error on the probability scale.
#'
#' @param doses Numeric vector of dose/threshold levels.
#' @param probs Numeric vector of true DLT probabilities at each dose level.
#' @param true_mtd Numeric. Threshold to treat as the true MTD (constraint point).
#' @param true_dlt_rate Numeric. DLT rate at the true MTD (constraint point).
#' @param link Character. Link function to use: "logit" or "probit".
#' @param plot Logical. Whether to plot the resulting curve.
#' @param plot_range Numeric vector of length 2. Range of thresholds to plot.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{true_prob_dlt}}{Function mapping thresholds to DLT probabilities.}
#'   \item{\code{alpha}}{Intercept of the model.}
#'   \item{\code{beta}}{Slope of the model.}
#'   \item{\code{fitted_probs}}{Model-predicted probabilities at input doses.}
#'   \item{\code{doses}}{Input dose levels.}
#'   \item{\code{probs}}{Input DLT probabilities.}
#'   \item{\code{true_mtd}}{The specified true MTD.}
#'   \item{\code{true_dlt_rate}}{The specified DLT rate at the true MTD.}
#'   \item{\code{link}}{Link function used.}
#' }
#'
#' @export
true_model <- function(
    doses,
    probs,
    true_mtd,
    true_dlt_rate,
    link = c("logit", "probit"),
    plot = TRUE,
    plot_range = range(doses)
) {
  stopifnot(length(doses) == length(probs))
  stopifnot(all(probs > 0 & probs < 1))
  
  link <- match.arg(link)
  qlink <- if (link == "logit") qlogis else qnorm
  plink <- if (link == "logit") plogis else pnorm
  
  # Link-scale value at the constraint
  target_link <- qlink(true_dlt_rate)
  
  # Define loss as function of beta only (alpha is solved from constraint)
  loss_fn <- function(beta) {
    alpha <- target_link - beta * true_mtd
    pred <- plink(alpha + beta * doses)
    sum((pred - probs)^2)
  }
  
  # Optimize beta under constraint
  opt <- optimize(loss_fn, interval = c(-100, 100))
  beta <- opt$minimum
  alpha <- target_link - beta * true_mtd
  
  # Final model function
  true_prob_dlt <- function(thresh) plink(alpha + beta * thresh)
  fitted_probs <- true_prob_dlt(doses)
  
  if (plot) {
    dose_seq <- seq(plot_range[1], plot_range[2], by = 0.001)
    plot(dose_seq, true_prob_dlt(dose_seq),
         type = "l", lwd = 2, col = "blue",
         ylab = "DLT Probability", xlab = "Threshold",
         main = paste("Constrained Dose-Toxicity Curve (", link, " link)", sep = ""))
    
    points(doses, probs, pch = 19, col = "red")
    
    abline(h = true_dlt_rate, v = true_mtd, col = "darkgreen", lty = 2)
    text(true_mtd, true_dlt_rate - 0.03,
         labels = paste0("True MTD = ", true_mtd,
                         "\nDLT Rate = ", true_dlt_rate),
         col = "darkgreen", pos = 4)
    
    legend("topleft",
           legend = c("Fitted Curve", "Specified Points", "True MTD Constraint"),
           col = c("blue", "red", "darkgreen"),
           lty = c(1, NA, 2),
           pch = c(NA, 19, NA),
           bty = "n")
  }
  
  list(
    true_prob_dlt = true_prob_dlt,
    alpha = alpha,
    beta = beta,
    fitted_probs = fitted_probs,
    doses = doses,
    probs = probs,
    true_mtd = true_mtd,
    true_dlt_rate = true_dlt_rate,
    link = link
  )
}

