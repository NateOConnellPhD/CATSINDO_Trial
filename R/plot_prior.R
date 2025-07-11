#' Plot Prior vs. True Dose–DLT Curve
#'
#' @param prior_model A list from `get_prior()` with fitted alpha, beta, priors, etc.
#' @param true_model A list from `true_model()` with `true_prob_dlt()` and optional `true_dlt_rate`, `true_mtd`
#' @param plot_range Range of doses to plot (x-axis)
#' @param target_dlt_rate Optional; horizontal line showing target toxicity
#' @param n_draws Number of draws for credible intervals
#'
#' @return A `ggplot2` object showing the prior mean curve with credible intervals
#'         and the true dose–DLT curve
#' @export
plot_prior <- function(
    prior_model,
    true_model,
    plot_range      = range(c(prior_model$doses, true_model$true_mtd)),
    target_dlt_rate = NULL,
    n_draws         = 1000
) {
  library(ggplot2)
  library(tibble)
  library(dplyr)
  
  # Inverse link
  linkinv <- if (prior_model$link == "logit") plogis else pnorm
  
  # Raw dose grid
  dose_seq <- seq(plot_range[1], plot_range[2], by = 0.001)
  
  # Scale to match prior_model's internal scaling
  sdose_seq <- (dose_seq - prior_model$dose_base) * prior_model$dose_scale_factor
  
  # Sample priors
  intercept_samps <- rnorm(n_draws, mean = prior_model$alpha, sd = prior_model$int_prior)
  slope_samps     <- rnorm(n_draws, mean = prior_model$beta,  sd = prior_model$slope_prior)
  
  pred_mat <- sapply(seq_len(n_draws), function(i) {
    linkinv(intercept_samps[i] + slope_samps[i] * sdose_seq)
  })
  
  prior_df <- tibble(
    dose  = dose_seq,
    mean  = rowMeans(pred_mat),
    lower = apply(pred_mat, 1, quantile, 0.025),
    upper = apply(pred_mat, 1, quantile, 0.975)
  )
  
  # True model values
  true_df <- tibble(
    dose     = dose_seq,
    dlt_prob = true_model$true_prob_dlt(dose_seq)
  )
  
  # Y-axis limit
  ymax     <- max(prior_df$upper, true_df$dlt_prob)
  ylim_max <- min(1, ymax + 0.05)
  
  # Target DLT line
  if (is.null(target_dlt_rate) && !is.null(prior_model$true_dlt_rate)) {
    target_dlt_rate <- prior_model$true_dlt_rate
  }
  
  # MTD line
  mtd_line <- prior_model$true_mtd %||% true_model$true_mtd
  
  # Plot
  p <- ggplot() +
    geom_ribbon(data = prior_df,
                aes(x = dose, ymin = lower, ymax = upper),
                fill = "darkred", alpha = 0.15) +
    geom_line(data = prior_df,
              aes(x = dose, y = mean, color = "Prior"),
              size = 1.2) +
    geom_line(data = true_df,
              aes(x = dose, y = dlt_prob, color = "True"),
              linetype = "dashed", size = 1.2)
  
  if (!is.null(mtd_line)) {
    p <- p + geom_vline(xintercept = mtd_line, linetype = "dotted", color = "gray40")
  }
  
  if (!is.null(target_dlt_rate)) {
    p <- p +
      geom_hline(yintercept = target_dlt_rate, linetype = "dotted", color = "gray40") +
      annotate("text",
               x = plot_range[1],
               y = target_dlt_rate + 0.03,
               label = paste0("Target DLT = ", target_dlt_rate),
               hjust = 0, vjust = 0,
               color = "gray40", size = 3.5)
  }
  
  p +
    scale_color_manual(name = "", values = c(Prior = "darkred", True = "blue")) +
    labs(
      title = paste0("Prior vs. True Dose–Response (", prior_model$link, ")"),
      x     = "Dose",
      y     = "DLT Probability"
    ) +
    scale_y_continuous(limits = c(0, ylim_max), breaks = scales::pretty_breaks(n = 5)) +
    theme_minimal()
}
