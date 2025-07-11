#' Run a Dose-Finding Trial Using Visit-Based Timing and Patient Heterogeneity
#'
#' Simulates and evaluates a dose-finding trial where dose escalation decisions are made at
#' fixed time intervals rather than after fixed-size cohorts. The trial uses longitudinal visit
#' data with staggered patient enrollment and accounts for patient-specific random effects
#' in the dose–toxicity relationship.
#'
#' @param visits Data frame of patient-visit records, such as from `sim_data()`. Must contain:
#'   \code{patient_id}, \code{visit_day}.
#' @param prior_info List from `get_prior()` with prior parameters and model configuration.
#' @param dlt_model List from `true_model()` containing the population model and patient-specific
#'   heterogeneity via \code{true_prob_dlt()} and \code{patient_effects()}.
#' @param target_dlt_rate Numeric. Target DLT probability used to define the MTD.
#' @param starting_dose Numeric. Initial threshold/dose value to assign for visits in the first interval.
#'   Must match one of the values in \code{prior_info$doses}.
#' @param decision_interval Integer. Time in days between MTD re-estimation updates. Default is 30.
#' @param min_data_per_update Integer. Minimum number of observed DLT outcomes required before
#'   updating the model. Prevents fitting with too little data. Default is 10.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{visits}}{Data frame of all visit-level records with DLT outcomes and assigned doses.}
#'   \item{\code{mtd_estimates}}{Tibble of MTD estimates over time (by decision interval).}
#'   \item{\code{final_model}}{Final fitted `brms` model with patient-level random effects.}
#'   \item{\code{dlt_model}}{The `true_model()` used for simulating DLT outcomes.}
#'   \item{\code{threshold_range}}{Fine grid of dose values for plotting.}
#'   \item{\code{dose_steps}}{The discrete dose levels considered for estimation.}
#'   \item{\code{dose_base}}{Value subtracted from doses so that min(dose) = 0.}
#'   \item{\code{true_mtd}}{The true MTD passed in from prior_info.}
#'   \item{\code{target_dlt_rate}}{The user-specified target DLT probability.}
#'   \item{\code{prior_info}}{List of prior settings used in the model.}
#'   \item{\code{link}}{The link function used ("logit" or "probit").}
#'   \item{\code{final_dose_response}}{Posterior means of DLT at each dose.}
#'   \item{\code{init_params}}{List of initial trial parameters for reproducibility.}
#' }
#'
#' @details
#' This function simulates how a dose-finding trial would proceed under a rolling enrollment
#' structure. It re-estimates the MTD at pre-defined intervals (e.g., every 30 days) using a
#' hierarchical Bayesian logistic or probit model fit with `brms`, including random intercepts
#' for patient-specific DLT risk.
#'
#' The function uses patient-specific latent effects from the \code{dlt_model$patient_effects()} function
#' to simulate heterogeneity in the DLT probabilities. At each interval, the estimated MTD is the
#' dose with posterior mean DLT probability closest to the \code{target_dlt_rate}.
#'
#' @seealso \code{\link{sim_data}}, \code{\link{get_prior}}, \code{\link{true_model}}
#' @export
run_trial <- function(
    visits,
    prior_info,
    dlt_model,
    target_dlt_rate,
    starting_dose         = NULL,
    decision_interval     = 30,
    min_data_per_update   = 10,
    safety_lead_in        = NULL,
    max_overdose_prob     = 0.25,
    overdose_eval_threshold = NULL
) {
  library(dplyr); library(brms); library(tibble)
  
  dose_steps      <- sort(unique(prior_info$doses))
  dose_base       <- prior_info$dose_base
  dose_scale      <- prior_info$dose_scale_factor
  min_raw         <- dose_base
  max_raw         <- max(dose_steps)
  
  if (is.null(starting_dose) || !(starting_dose %in% dose_steps)) {
    stop("starting_dose must be one of prior_info$doses")
  }
  
  if (is.null(overdose_eval_threshold)) {
    overdose_eval_threshold <- target_dlt_rate
  }
  
  visits <- as.data.frame(ungroup(visits))
  visits$assigned_threshold <- NA_real_
  visits$realized_threshold <- NA_real_
  visits$dlt                 <- NA_integer_
  
  if (!"patient_id" %in% names(visits)) {
    stop("visits must include a 'patient_id' column")
  }
  
  patient_ids <- unique(visits$patient_id)
  patient_offsets <- dlt_model$patient_effects(length(patient_ids))
  names(patient_offsets) <- patient_ids
  
  lead_in_interval <- if (is.null(safety_lead_in)) decision_interval else safety_lead_in
  
  first_batch <- visits$visit_day <= lead_in_interval
  visits$realized_threshold[first_batch] <- starting_dose
  visits$assigned_threshold[first_batch] <- (starting_dose - dose_base) * dose_scale
  visits$dlt[first_batch] <- mapply(
    function(pid) rbinom(1, 1, dlt_model$true_prob_dlt(starting_dose, pid, patient_offsets)),
    visits$patient_id[first_batch]
  )
  
  current_time <- lead_in_interval
  mtd_estimates <- tibble(visit_day = 0, mtd = starting_dose)
  
  priors <- c(
    set_prior(
      paste0("normal(", prior_info$beta,  ",", prior_info$slope_prior, ")"),
      class = "b", lb = 0
    ),
    set_prior(
      paste0("normal(", prior_info$alpha, ",", prior_info$int_prior, ")"),
      class = "Intercept"
    ),
    set_prior("normal(0, 0.3)", class = "sd", group = "patient_id", lb = 0)
  )
  
  fit <- NULL
  max_time <- max(visits$visit_day)
  
  while (current_time <= max_time) {
    df <- filter(visits, visit_day <= current_time, !is.na(dlt))
    
    if (nrow(df) < min_data_per_update) {
      current_time <- current_time + decision_interval
      next
    }
    
    fit <- brm(
      dlt ~ assigned_threshold + (1 | patient_id),
      data    = df,
      family  = bernoulli(link = prior_info$link),
      prior   = priors,
      chains  = 4, iter = 2000,
      refresh = 0, silent = TRUE,
      backend = "cmdstanr"
    )
    
    step_df <- data.frame(assigned_threshold = (dose_steps - dose_base) * dose_scale)
    step_post_mat <- posterior_epred(fit, newdata = step_df, re_formula = NA)
    step_mean <- colMeans(step_post_mat)
    step_over_prob <- colMeans(step_post_mat > overdose_eval_threshold)
    
    safe_doses <- which(step_over_prob <= max_overdose_prob)
    
    if (length(safe_doses) > 0) {
      idx_best <- safe_doses[which.min(abs(step_mean[safe_doses] - target_dlt_rate))]
    } else {
      idx_best <- 1
    }
    
    nominal_mtd <- dose_steps[idx_best]
    
    idx_next <- which(
      visits$visit_day > current_time &
        visits$visit_day <= (current_time + decision_interval) &
        is.na(visits$realized_threshold)
    )
    
    visits$realized_threshold[idx_next] <- nominal_mtd
    visits$assigned_threshold[idx_next] <- (nominal_mtd - dose_base) * dose_scale
    visits$dlt[idx_next] <- mapply(
      function(pid) rbinom(1, 1, dlt_model$true_prob_dlt(nominal_mtd, pid, patient_offsets)),
      visits$patient_id[idx_next]
    )
    
    mtd_estimates <- bind_rows(mtd_estimates, tibble(visit_day = current_time, mtd = nominal_mtd))
    current_time <- current_time + decision_interval
  }
  
  step_pred_df <- data.frame(assigned_threshold = (dose_steps - dose_base) * dose_scale)
  step_post_mat <- posterior_epred(fit, newdata = step_pred_df, re_formula = NA)
  
  final_dose_response <- tibble(
    threshold = dose_steps,
    mean_dlt  = colMeans(step_post_mat),
    lower_ci  = apply(step_post_mat, 2, quantile, 0.025),
    upper_ci  = apply(step_post_mat, 2, quantile, 0.975),
    over_prob = colMeans(step_post_mat > overdose_eval_threshold)
  )
  
  final_mtd <- tail(mtd_estimates$mtd, 1)
  
  list(
    visits               = visits,
    mtd_estimates        = mtd_estimates,
    final_mtd            = final_mtd,  
    final_model          = fit,
    dlt_model            = dlt_model,
    dose_steps           = dose_steps,
    true_mtd             = dlt_model$true_mtd,
    target_dlt_rate      = target_dlt_rate,
    prior_info           = prior_info,
    link                 = prior_info$link,
    final_dose_response  = final_dose_response,
    init_params = list(
      decision_interval       = decision_interval,
      safety_lead_in          = safety_lead_in,
      starting_dose           = starting_dose,
      min_data_per_update     = min_data_per_update,
      max_overdose_prob       = max_overdose_prob,
      overdose_eval_threshold = overdose_eval_threshold
    )
  )
}

set.seed(2025)
result <- run_noninf_sim(n_patients = 40, true_dlt_rate = 0.15, ni_margin = 0.20)

# Plot
library(ggplot2)
ggplot(result, aes(x = day, y = power)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  labs(title = "Power to Declare Non-Inferiority Over Time",
       x = "Day",
       y = "Estimated Power") +
  theme_minimal()

