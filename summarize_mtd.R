#’ Summarize Final Estimated MTD Across Simulations
#’
#’ Given a list of results from \code{run_simulations()}, extract each simulation’s
#’ final MTD estimate and compute overall summary statistics (mean, median, SD, IQR)
#’ plus a frequency distribution of the selected doses.
#’
#’ @param results_list A list of length \code{n_sim}, each element a result list from \code{run_simulations()}.
#’ @param mtd_field Character. The name of the element within \code{res$mtd_estimates} that holds the numeric MTD estimate. Defaults to \code{"estimate"}.
#’ @return A list with two elements:
#’   \describe{
#’     \item{stats}{Named numeric vector: mean, median, sd, 25th and 75th percentiles of the estimated MTDs.}
#’     \item{freq}{Table of how often each dose level was chosen as the final MTD.}
#’   }
#’ @export
summarize_mtd <- function(results_list, mtd_field = "estimate") {
  
  # Extract the MTD estimate from each simulation
  mtd_vec <- vapply(
    results_list,
    function(res) {
      if (!is.list(res$mtd_estimates) || is.null(res$mtd_estimates[[mtd_field]])) {
        stop("Cannot find field '", mtd_field, "' in one of the result elements.")
      }
      res$mtd_estimates[[mtd_field]]
    },
    numeric(1)
  )
  
  # Compute summary statistics
  stats <- c(
    mean   = mean(mtd_vec),
    median = median(mtd_vec),
    sd     = sd(mtd_vec),
    q25    = as.numeric(quantile(mtd_vec, 0.25)),
    q75    = as.numeric(quantile(mtd_vec, 0.75))
  )
  
  # Frequency table of chosen MTDs
  freq <- table(mtd_vec)
  
  list(
    stats = stats,
    freq  = freq
  )
}
