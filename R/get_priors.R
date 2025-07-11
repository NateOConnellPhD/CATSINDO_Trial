#' Derive Priors for Slope and Intercept for Logistic/Probit Model Parameters
#'
#' Given two points on the dose-toxicity curve, this function calculates the implied
#' prior means for the intercept and slope of a logistic or probit regression model,
#' where the predictor is a scaled dose (to [0, 1]). Also accepts and returns SDs for priors.
#'
#' @param anchor_thresh Numeric. Dose where DLT probability is known (e.g., 0.10).
#' @param anchor_dlt_rate Numeric. DLT probability at \code{anchor_thresh} (e.g., 0.10).
#' @param est_mtd Numeric. Dose assumed to be the MTD (e.g., 0.18).
#' @param target_dlt_rate Numeric. Target DLT probability at the estimated MTD (e.g., 0.15).
#' @param thresh_range Numeric vector of length 2. The minimum and maximum dose values for scaling.
#' @param prior_sds Named list with elements `intercept` and `slope` specifying prior SDs.
#' @param link Character. Link function to use: "logit" (default) or "probit".
#'
#' @return A named list with:
#' \describe{
#'   \item{\code{intercept_mean}}{Prior mean for the intercept (alpha).}
#'   \item{\code{slope_mean}}{Prior mean for the slope (beta).}
#'   \item{\code{link}}{Link function used: "logit" or "probit".}
#'   \item{\code{x_anchor}}{Scaled anchor threshold (between 0 and 1).}
#'   \item{\code{x_target}}{Scaled estimated MTD (between 0 and 1).}
#'   \item{\code{anchor}}{Unscaled anchor threshold.}
#'   \item{\code{target}}{Unscaled estimated MTD.}
#'   \item{\code{thresh_range}}{Original threshold range used for scaling.}
#'   \item{\code{prior_sds}}{List containing `intercept` and `slope` SDs.}
#' }
#'
#' @examples
#' get_priors(
#'   anchor_thresh = 0.10, anchor_dlt_rate = 0.10,
#'   est_mtd = 0.18, target_dlt_rate = 0.15,
#'   thresh_range = c(0.10, 0.25),
#'   prior_sds = list(intercept = 0.5, slope = 2)
#' )
#'
#' @export
get_priors_old <- function(
    anchor_thresh, 
    anchor_dlt_rate,
    est_mtd, 
    target_dlt_rate,
    thresh_range, 
    link = c("logit", "probit"),
    prior_sds = list(intercept = 0.5, slope = 2)
) {
  link <- match.arg(link)
  qlink <- if (link == "logit") qlogis else qnorm
  
  # 1) compute scaled-x for each point
  dose_range <- diff(thresh_range)           # T₂ – T₁
  x_anchor  <- (anchor_thresh - thresh_range[1]) / dose_range
  x_target  <- (est_mtd      - thresh_range[1]) / dose_range
  
  # 2) map your two probs onto the link scale
  link_anchor <- qlink(anchor_dlt_rate)
  link_target <- qlink(target_dlt_rate)
  
  # 3) get the slope *per unit scaled‑x*
  slope_scaled <- (link_target - link_anchor) / (x_target - x_anchor)
  
  # 4) convert back to “per unit raw dose”
  slope_raw     <- slope_scaled / dose_range
  intercept_raw <- link_anchor - slope_raw * anchor_thresh
  
  list(
    intercept_mean  = as.numeric(intercept_raw),
    slope_mean      = as.numeric(slope_raw),
    intercept_sd    = prior_sds$intercept,
    slope_sd        = prior_sds$slope,
    link            = link,
    anchor_thresh   = anchor_thresh,
    anchor_dlt_rate = anchor_dlt_rate,
    est_mtd         = est_mtd,
    target_dlt_rate = target_dlt_rate,
    thresh_range    = thresh_range
  )
}
