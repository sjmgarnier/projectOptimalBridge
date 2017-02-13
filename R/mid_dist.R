#' @title Mid-Distance Function
#'
#' @description This function computes the distance between two bridges along a
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
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @export
#'
mid_dist <- function(d1, d2, theta, Dmax) {
  args <- as.list(environment())

  if (!all(sapply(args, is.numeric)))
    stop("All arguments must be numeric.")

  l <- lengths(args)
  if (!all(((l == 1) + (l == max(l))) > 0))
    stop("All arguments should be of either length 1 or of the same length.")

  abs((1 / cos(0.5 * theta)) * (Dmax - d1 - d2))
}
