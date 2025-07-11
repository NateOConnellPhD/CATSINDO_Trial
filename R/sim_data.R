#' Simulate Longitudinal HDMTX Visits for Staggered Patient Enrollment
#'
#' Each patient is enrolled according to a Poisson process (exponential gaps),
#' and receives a fixed 12-visit HDMTX schedule based on their personal enrollment day.
#'
#' @param n_patients Integer. Number of patients to simulate.
#' @param mean_enrollment_gap Numeric. Mean days between enrollments (Exponential distribution).
#'
#' @return A data frame with one row per patient-visit combination:
#' \describe{
#'   \item{patient_id}{Patient ID.}
#'   \item{enrollment_day}{Day patient enrolled (study day).}
#'   \item{visit_number}{1 through 12.}
#'   \item{visit_day}{Absolute day of this visit since study start.}
#'   \item{days_since_enrollment}{Relative day of visit for the patient.}
#' }
#'
#' @examples
#' set.seed(123)
#' sim_data(n_patients = 3, mean_enrollment_gap = 7)
#'
#' @export
# sim_data <- function(n_patients,
#                      mean_enrollment_gap = 7) {
#   library(dplyr)
#   
#   # HDMTX visit schedule in days relative to each patient’s enrollment
#   mtx_weeks <- c(1, 2, 6, 7, 11, 12, 16, 17, 20, 21, 24, 25)
#   visit_offsets <- (mtx_weeks - 1) * 7 + 1  # e.g., Week 1 = Day 1, Week 2 = Day 8, etc.
#   
#   # Simulate enrollment times using exponential gaps
#   enrollment_gaps <- rexp(n_patients, rate = 1 / mean_enrollment_gap)
#   enrollment_days <- round(cumsum(enrollment_gaps))
#   
#   # Expand to patient-visits
#   visits <- lapply(1:n_patients, function(pid) {
#     enroll_day <- enrollment_days[pid]
#     tibble(
#       patient_id = pid,
#       enrollment_day = enroll_day,
#       visit_number = 1:12,
#       days_since_enrollment = visit_offsets,
#       visit_day = enroll_day + visit_offsets
#     )
#   }) %>%
#     bind_rows() %>%
#     arrange(visit_day)
#   
#   return(visits)
# }

#' @examples
#' set.seed(123)
#' sim_data(n_patients = 3, mean_enrollment_gap = 7)
#'
#' @export
sim_data <- function(n_patients,
                     mean_enrollment_gap = 7,
                     mean_visits = 10,
                     max_visits = 12) {
  library(dplyr)
  
  # HDMTX visit schedule in days relative to enrollment
  mtx_weeks <- c(1, 2, 6, 7, 11, 12, 16, 17, 20, 21, 24, 25)
  visit_offsets_all <- (mtx_weeks - 1) * 7 + 1
  
  # Simulate enrollment times using exponential gaps
  enrollment_gaps <- rexp(n_patients, rate = 1 / mean_enrollment_gap)
  enrollment_days <- round(cumsum(enrollment_gaps))
  
  # Simulate number of visits per patient with a cap at max_visits
  num_visits <- rpois(n_patients, lambda = mean_visits)
  num_visits <- pmin(num_visits, max_visits)  # cap at 12
  num_visits[num_visits < 1] <- 1            # ensure at least 1 visit
  
  # Force ~15–20% to have ≤6 visits
  low_visit_ids <- sample(1:n_patients, size = round(0.17 * n_patients))
  num_visits[low_visit_ids] <- sample(3:6, length(low_visit_ids), replace = TRUE)
  
  # Expand to patient-visits
  visits <- lapply(1:n_patients, function(pid) {
    n_visits <- num_visits[pid]
    offsets <- visit_offsets_all[1:n_visits]
    enroll_day <- enrollment_days[pid]
    
    tibble(
      patient_id = pid,
      enrollment_day = enroll_day,
      visit_number = 1:n_visits,
      days_since_enrollment = offsets,
      visit_day = enroll_day + offsets
    )
  }) %>%
    bind_rows() %>%
    arrange(visit_day)
  
  return(visits)
}

