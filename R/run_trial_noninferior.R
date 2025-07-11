#' Run Dose-Finding Trial Using Non-Inferiority Based Escalation
#'
#' Selects doses where the 95% credible interval includes the benchmark DLT rate,
#' and escalates to the highest such dose. The final MTD is the highest dose
#' whose posterior CrI contains the benchmark rate (e.g., 0.10).
#'
#' @param visits Data frame of patient visit data.
#' @param prior_info Output from `get_prior()`, including priors and dose scale.
#' @param dlt_model Output from `true_model()` with true_prob_dlt().
#' @param target_dlt_rate Benchmark DLT rate (e.g., 0.10).
#' @param cohort_size Number of visits per cohort after the first.
#' @param first_cohort_size Number of visits in the initial cohort.
#' @param starting_dose Raw dose value for first cohort (must match prior_info$doses).
#'
#' @return A list with trial results, including visit data, estimated dose-response curve,
#'         and the final MTD.
#' @export
run_trial_noninferiority <- function(
    visits,
    prior_info,
    dlt_model,
    target_dlt_rate,
    cohort_size = 10,
    first_cohort_size = 10,
    starting_dose = NULL
) {
  library(dplyr); library(brms); library(tibble)
  
  dose_steps      <- sort(unique(prior_info$doses))
  dose_base       <- prior_info$dose_base
  threshold_range <- seq(dose_base, max(dose_steps), by = 0.001)
  
  if (is.null(starting_dose) || !(starting_dose %in% dose_steps)) {
    stop("starting_dose must be one of prior_info$doses")
  }
  
  visits <- as.data.frame(ungroup(visits))
  visits$assigned_threshold <- NA_real_
  visits$realized_threshold <- NA_real_
  visits$dlt <- NA_integer_
  
  visits$assigned_threshold[1:first_cohort_size] <- starting_dose - dose_base
  visits$realized_threshold[1:first_cohort_size] <- starting_dose
  visits$dlt[1:first_cohort_size] <- rbinom(
    first_cohort_size, 1,
    dlt_model$true_prob_dlt(starting_dose)
  )
  
  mtd_estimates <- tibble(visit = integer(), mtd = numeric())
  
  priors <- c(
    set_prior(
      paste0("normal(", prior_info$beta,  ",", prior_info$slope_prior, ")"),
      class = "b", lb = 0
    ),
    set_prior(
      paste0("normal(", prior_info$alpha, ",", prior_info$int_prior, ")"),
      class = "Intercept"
    ),
    set_prior("exponential(1)", class = "sd", group = "patient_id")
  )
  
  fit <- NULL
  for (i in seq(first_cohort_size, nrow(visits), by = cohort_size)) {
    df <- visits[1:i, ] %>% filter(!is.na(dlt))
    
    fit <- brm(
      dlt ~ assigned_threshold + (1 | patient_id),
      data    = df,
      family  = bernoulli(link = prior_info$link),
      prior   = priors,
      chains  = 4, iter = 2000,
      refresh = 0, silent = TRUE,
      backend = "cmdstanr"
    )
    
    # Posterior prediction at discrete dose levels
    step_df <- data.frame(assigned_threshold = dose_steps - dose_base)
    step_mat <- posterior_epred(fit, newdata = step_df, re_formula = NA)
    step_mean <- colMeans(step_mat)
    step_ci <- apply(step_mat, 2, quantile, probs = c(0.025, 0.975))
    
    # Find the highest dose where CrI includes the target rate
    eligible_doses <- which(step_ci[1, ] <= target_dlt_rate & step_ci[2, ] >= target_dlt_rate)
    
    if (length(eligible_doses) > 0) {
      nominal_mtd <- dose_steps[max(eligible_doses)]
    } else {
      nominal_mtd <- starting_dose
    }
    
    mtd_estimates <- bind_rows(mtd_estimates, tibble(visit = i, mtd = nominal_mtd))
    
    idx_next <- (i + 1):min(i + cohort_size, nrow(visits))
    visits$realized_threshold[idx_next] <- nominal_mtd
    visits$assigned_threshold[idx_next] <- nominal_mtd - dose_base
    visits$dlt[idx_next] <- rbinom(length(idx_next), 1, dlt_model$true_prob_dlt(nominal_mtd))
  }
  
  # Final dose-response summaries
  step_post_final <- posterior_epred(fit, newdata = step_df, re_formula = NA)
  final_dose_response <- tibble(
    threshold = dose_steps,
    mean_dlt = colMeans(step_post_final),
    lower_95 = apply(step_post_final, 2, quantile, 0.025),
    upper_95 = apply(step_post_final, 2, quantile, 0.975)
  )
  
  list(
    visits = visits,
    mtd_estimates = mtd_estimates,
    final_model = fit,
    dlt_model = dlt_model,
    threshold_range = threshold_range,
    dose_steps = dose_steps,
    dose_base = dose_base,
    true_mtd = prior_info$true_mtd,
    target_dlt_rate = target_dlt_rate,
    prior_info = prior_info,
    link = prior_info$link,
    final_dose_response = final_dose_response,
    init_params = list(
      first_cohort_size = first_cohort_size,
      cohort_size = cohort_size,
      starting_dose = starting_dose
    )
  )
}
