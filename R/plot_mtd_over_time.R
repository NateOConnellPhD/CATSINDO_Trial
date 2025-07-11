#' Plot Realized MTD Over Time with Posterior DLT Estimates
#'
#' Plots the assigned threshold (MTD) across visit days from a run_trial() output,
#' and adds posterior mean DLT estimates with 95% credible intervals on the right axis.
#'
#' @param results A list returned by `run_trial()`
#'
#' @return A ggplot object
#' @export
plot_mtd_over_time <- function(results) {
  library(ggplot2)
  library(dplyr)
  
  final_mtd <- tail(results$mtd_estimates$mtd, 1)
  x_max <- max(results$visits$visit_day)
  
  ggplot(results$visits, aes(x = visit_day, y = realized_threshold)) +
    geom_point(alpha = 0.3, size = 1) +
    geom_step(data = results$mtd_estimates, aes(x = visit_day, y = mtd),
              color = "blue", linewidth = 1.2, direction = "hv") +
    geom_hline(yintercept = results$true_mtd, linetype = "dashed", color = "red") +
    geom_hline(yintercept = final_mtd, linetype = "dotted", color = "blue") +
    geom_text(
      data = results$final_dose_response,
      aes(
        x = x_max + 40,
        y = threshold,
        label = sprintf("%.2f (%.2f, %.2f)", mean_dlt, lower_ci, upper_ci)
      ),
      hjust = 0,
      nudge_x = 130,
      size = 3.3
    ) +
    coord_cartesian(xlim = c(0, x_max + 120), clip = "off") +
    scale_y_continuous(
      results$dose_steps,
      name = "Dose Steps"
    ) +
    labs(
      title = "Estimated MTD Over Time with Posterior DLT Rates (Right Side)",
      x = "Visit Day"
    ) +
    theme_minimal() +
    theme(
      plot.margin = margin(t = 10, r = 140, b = 10, l = 10)
    )
}
