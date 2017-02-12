#' @title Mid-Distance Function
#'
#' @description This function models the distance between two bridges along a
#'  common arm.
#'
#' @param d1 A numeric vector of distances of first bridge from apex in cm.
#'
#' @param d2 A numeric vector of distances of second bridge from apex in cm.
#'
#' @param theta A numeric vector of angles of apex in rad.
#'
#' @param Dmax A numeric vector of the maximum possible distances of a bridge
#'  from apex in cm.
#'
#' @return A vector of the same length as \code{d1, d2, theta}.
#'
#' @author Jason Graham, \email{jason.graham@@scranton.edu}
#'
#' @export
#'
mid_dist <- function(d1, d2, theta, Dmax) {
  if (!all(is.numeric(d1), is.numeric(d2), is.numeric(theta), is.numeric(Dmax)))
    stop("All arguments must be numeric")

  l <- lengths(list(d1, d2, theta, Dmax))
  if (!all(((l == 1) + (l == max(l))) > 0))
    stop("All arguments should be of either length 1 or of the same length.")

  abs((1 / cos(0.5 * theta)) * (Dmax - d1 - d2))
}
