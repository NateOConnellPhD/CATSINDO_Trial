#' Plot Estimated, True, and Prior Dose-Response Curves
#'
#' @param results Output list from `run_trial()`
#' @return A ggplot object
#' @export
plot_est_dr <- function(results) {
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(brms)
  
  # Extract
  final_model     <- results$final_model
  mtd_estimates   <- results$mtd_estimates
  threshold_range <- results$threshold_range
  prior_info      <- results$prior_info
  dlt_model       <- results$dlt_model
  dose_base       <- results$dose_base
  dose_scale      <- prior_info$dose_scale_factor
  
  min_thresh <- min(threshold_range)
  max_thresh <- max(threshold_range)
  
  # 1) Estimated (posterior) dose-response
  pred_df <- data.frame(
    assigned_threshold = (threshold_range - dose_base) * dose_scale
  )
  post_mat <- posterior_epred(final_model,
                              newdata    = pred_df,
                              re_formula = NA)
  
  est_df <- tibble(
    threshold = threshold_range,
    mean      = colMeans(post_mat),
    lower     = apply(post_mat, 2, quantile, probs = 0.025),
    upper     = apply(post_mat, 2, quantile, probs = 0.975)
  )
  
  # 2) True dose-response
  true_df <- tibble(
    threshold = threshold_range,
    true_prob = dlt_model$true_prob_dlt(threshold_range)
  )
  
  # 3) Prior dose-response curve
  link_fn <- if (prior_info$link == "probit") pnorm else plogis
  prior_probs <- link_fn(
    prior_info$alpha + prior_info$beta * ((threshold_range - dose_base) * dose_scale)
  )
  prior_df <- tibble(
    threshold = threshold_range,
    y         = prior_probs,
    curve     = "Prior"
  )
  
  # 4) Combine all three
  plot_curves <- bind_rows(
    est_df  %>% select(threshold, y = mean) %>% mutate(curve = "Estimated"),
    true_df %>% rename(y = true_prob)       %>% mutate(curve = "True"),
    prior_df
  )
  
  # 5) Vertical MTD lines
  final_mtd <- as.numeric(tail(mtd_estimates$mtd, 1))
  true_mtd  <- results$true_mtd
  
  # 6) Target DLT threshold
  if (!is.null(results$init_params$target_type) &&
      results$init_params$target_type == "relative") {
    target_prob <- dlt_model$anchor_dlt_rate + results$target_dlt_rate
  } else {
    target_prob <- results$target_dlt_rate
  }
  
  # 7) Plot
  ggplot() +
    geom_ribbon(data = est_df,
                aes(x = threshold, ymin = lower, ymax = upper),
                fill = "steelblue", alpha = 0.2) +
    geom_line(data = plot_curves,
              aes(x = threshold, y = y, color = curve, linetype = curve),
              size = 1.2) +
    geom_hline(yintercept = target_prob,
               linetype = "dotted", color = "red", size = 1) +
    annotate("text",
             x     = min_thresh,
             y     = target_prob + 0.02,
             label = paste0("Target DLT = ", round(target_prob, 3)),
             hjust = 0, color = "red") +
    geom_vline(xintercept = final_mtd,
               linetype = "dotdash", color = "purple", size = 1) +
    geom_vline(xintercept = true_mtd,
               linetype = "dashed",  color = "black",  size = 1) +
    annotate("text",
             x     = final_mtd,
             y     = 0.05,
             label = paste0("Final MTD = ", round(final_mtd, 3)),
             angle = 90, vjust = -0.5, color = "purple") +
    annotate("text",
             x     = true_mtd,
             y     = 0.05,
             label = paste0("True MTD = ",  round(true_mtd,  3)),
             angle = 90, vjust =  1.5, color = "black") +
    labs(
      title    = "Estimated, True, and Prior Dose‑Response Curves",
      x        = "Threshold (raw dose)",
      y        = "Probability of DLT",
      color    = "Curve",
      linetype = "Curve"
    ) +
    scale_color_manual(values = c(
      "Estimated" = "steelblue",
      "True"      = "darkgreen",
      "Prior"     = "darkorange"
    )) +
    scale_linetype_manual(values = c(
      "Estimated" = "solid",
      "True"      = "dashed",
      "Prior"     = "dotdash"
    )) +
    scale_x_continuous(limits = c(min_thresh, max_thresh)) +
    theme_minimal(base_size = 14)
}

