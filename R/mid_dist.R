#' @title Mid-distance function
#'
#' @description This function models the distance between two
#'              bridges along a common arm.
#'
#' @param d1 Distance of one bridge from apex.
#'
#' @param d2 Distance of other bridge from apex.
#'
#' @param theta Angle of apex (in radians)
#'
#' @param Dmax Maximum possible distance of a bridge from apex.
#'
#' @return mid-distance
#'
#' @author Jason Graham, \email{jason.graham@@scranton.edu}
#'
#' @export
#'
mid_dist <- function(d1, d2, theta, Dmax){

  return(abs((1/cos(0.5*theta))*(Dmax - d1 - d2)))

}
