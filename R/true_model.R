#' Construct a true dose-toxicity model with optional patient-level variation
#'
#' @param doses Vector of dose levels
#' @param probs Vector of DLT probabilities at each dose
#' @param target_dlt_rate DLT probability to enforce at the MTD (default 0.15)
#' @param plot Whether to generate a plot (default FALSE)
#' @param plot_range Optional range for plotting
#' @param sigma_re SD of random intercepts for patient-level effects (default 1)
#' @param n_patients_spaghetti Number of spaghetti lines to simulate (default 50)
#' @param dose_scale_factor Multiplicative scaling factor applied after shifting doses to min(doses)
#'
#' @return A list with alpha, beta, dose_base, dose_scale_factor, patient_effects,
#'         true_prob_dlt, true_mtd, and true_dlt_rate
#'
#' @export
true_model <- function(doses, probs, target_dlt_rate = 0.15,
                       plot = FALSE, plot_range = NULL, sigma_re = 1,
                       n_patients_spaghetti = 50, dose_scale_factor = 10) {
  library(ggplot2)
  library(dplyr)
  
  # Scale doses
  dose_base <- min(doses)
  sdoses <- (doses - dose_base) * dose_scale_factor
  
  # Find MTD index and scaled MTD
  idx_mtd <- which.min(abs(probs - target_dlt_rate))
  smtd <- sdoses[idx_mtd]
  logit_target <- qlogis(target_dlt_rate)
  
  # Optimize slope with intercept fixed to match MTD
  loss_fn <- function(beta) {
    alpha <- logit_target - beta * smtd
    logits <- qlogis(probs)
    sum((logits - (alpha + beta * sdoses))[-idx_mtd]^2)
  }
  
  beta <- optimize(loss_fn, interval = c(-10, 10))$minimum
  alpha <- logit_target - beta * smtd
  
  # Population DLT curve
  population_dlt <- function(dose) {
    scaled_dose <- (dose - dose_base) * dose_scale_factor
    plogis(alpha + beta * scaled_dose)
  }
  
  # Recalculate MTD and DLT rate using final curve
  pred_dlt_probs <- population_dlt(doses)
  idx_mtd <- which.min(abs(pred_dlt_probs - target_dlt_rate))
  true_mtd <- doses[idx_mtd]
  true_dlt_rate <- pred_dlt_probs[idx_mtd]
  
  # Patient-specific intercepts
  patient_effects <- function(n_patients = n_patients_spaghetti) {
    rnorm(n_patients, mean = 0, sd = sigma_re)
  }
  
  # Dose–DLT with optional patient-level offset
  true_prob_dlt <- function(dose, patient_id = NULL, patient_offsets = NULL) {
    scaled_dose <- (dose - dose_base) * dose_scale_factor
    lin_pred <- alpha + beta * scaled_dose
    if (!is.null(patient_id) && !is.null(patient_offsets)) {
      lin_pred <- lin_pred + patient_offsets[patient_id]
    }
    plogis(lin_pred)
  }
  
  # Plot if requested
  if (plot) {
    plot_seq <- if (is.null(plot_range)) seq(min(doses), max(doses), by = 0.001) else seq(plot_range[1], plot_range[2], by = 0.001)
    preds <- population_dlt(plot_seq)
    
    offsets <- patient_effects(n_patients_spaghetti)
    spaghetti_data <- lapply(1:n_patients_spaghetti, function(id) {
      data.frame(
        dose = plot_seq,
        dlt_prob = true_prob_dlt(plot_seq, patient_id = rep(id, length(plot_seq)), patient_offsets = offsets),
        patient_id = id
      )
    }) %>% bind_rows()
    
    p <- ggplot() +
      geom_line(data = spaghetti_data, aes(x = dose, y = dlt_prob, group = patient_id),
                color = "gray70", alpha = 0.5) +
      geom_point(aes(x = doses, y = probs), color = "red", size = 3) +
      geom_line(aes(x = plot_seq, y = preds), color = "blue", size = 1) +
      geom_hline(yintercept = true_dlt_rate, linetype = "dashed", color = "darkred") +
      geom_vline(xintercept = true_mtd, linetype = "dotted", color = "darkred") +
      labs(
        x = "Dose", y = "True DLT Probability",
        title = "True Dose–DLT Curve with Population and Patient-Level Effects"
      ) +
      theme_minimal()
    
    print(p)
  }
  
  list(
    alpha = alpha,
    beta = beta,
    dose_base = dose_base,
    dose_scale_factor = dose_scale_factor,
    patient_effects = patient_effects,
    true_prob_dlt = true_prob_dlt,
    true_mtd = true_mtd,
    true_dlt_rate = true_dlt_rate
  )
}


# 
# true_model <- function(
#     doses,
#     probs,
#     true_mtd = NULL,
#     true_dlt_rate = NULL,
#     plot = TRUE,
#     plot_range = range(doses)
# ) {
#   stopifnot(length(doses) == length(probs))
#   stopifnot(all(probs > 0 & probs < 1))
#   
#   # Spline interpolation function
#   smooth_fn <- splinefun(doses, probs, method = "monoH.FC")
#   
#   if (plot) {
#     dose_seq <- seq(plot_range[1], plot_range[2], by = 0.001)
#     plot(dose_seq, smooth_fn(dose_seq),
#          type = "l", lwd = 2, col = "blue",
#          ylab = "DLT Probability", xlab = "Threshold",
#          main = "True Dose-Toxicity Curve")
#     
#     points(doses, probs, pch = 19, col = "red")
#     
#     if (!is.null(true_mtd) && !is.null(true_dlt_rate)) {
#       abline(h = true_dlt_rate, v = true_mtd, col = "darkgreen", lty = 2)
#       text(true_mtd, true_dlt_rate - 0.03,
#            labels = paste0("True MTD = ", true_mtd,
#                            "\nDLT Rate = ", true_dlt_rate),
#            col = "darkgreen", pos = 4)
#     }
#     
#     # legend("topleft",
#     #        legend = c("Smooth True Curve", "Defined Points",
#     #                   if (!is.null(true_mtd) && !is.null(true_dlt_rate)) "True MTD Constraint" else NULL),
#     #        col = c("blue", "red", "darkgreen")[1:(2 + !is.null(true_mtd) && !is.null(true_dlt_rate))],
#     #        lty = c(1, NA, 2)[1:(2 + !is.null(true_mtd) && !is.null(true_dlt_rate))],
#     #        pch = c(NA, 19, NA)[1:(2 + !is.null(true_mtd) && !is.null(true_dlt_rate))],
#     #        bty = "n")
#   }
#   
#   list(
#     true_prob_dlt = smooth_fn,
#     doses = doses,
#     probs = probs,
#     true_mtd = true_mtd,
#     true_dlt_rate = true_dlt_rate
#   )
# }
