#' Scale values to [0, 1] based on custom min and max
#'
#' @param x Numeric vector to scale.
#' @param min_val Numeric. The minimum value of the original scale.
#' @param max_val Numeric. The maximum value of the original scale.
#'
#' @return A numeric vector scaled to [0, 1].
#' @export
scale_dose <- function(x, min_val, max_val) {
  (x - min_val) / (max_val - min_val)
}


#' Inverse scale from [0, 1] back to original scale
#'
#' @param x Numeric vector in [0, 1] to unscale.
#' @param min_val Numeric. The minimum value of the original scale.
#' @param max_val Numeric. The maximum value of the original scale.
#'
#' @return A numeric vector mapped back to the original scale.
#' @export
inv_scale_dose  <- function(x, min_val, max_val) {
  x * (max_val - min_val) + min_val
}
