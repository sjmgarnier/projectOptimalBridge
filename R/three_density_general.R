#' @title Density Function for Three-Apex Apparatus
#'
#' @description This function computes the density as a function of distance for
#'  the three-apex appratus.
#'
#' @param d A length 3 numeric dector of distances values in cm.
#'
#' @param angle A scalar for the angle of apex in deg.
#'
#' @param ant.dens A scalar for the density of ants in ants/cm (default: 2.2).
#'
#' @param LT A scalar for the overall trail length in cm (default: 100).
#'
#' @param L0 A scalar for the length of the apparatus arm in cm (default: 22).
#'
#' @param wA A scalar for the width of the arms of the apparatus in cm (default: 3.3).
#'
#' @param ln A scalar for the average length of an ant in cm (default:0.691).
#'
#' @param wn A scalar for the average width of an ant in cm (default: 0.107).
#'
#' @param alpha A scalar for the free fitting parameter value from Reid et al.
#'  (default: 17.02).
#'
#' @return A vector of the of the same length as \code{d, angle}.
#'
#' @author Jason Graham, \email{jason.graham@@scranton.edu}
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @export
#'
three_density_general <- function(d, angle, ant.dens = 2.2, LT = 100, L0 = 22,
                                  wA = 3.3, ln = 0.691, wn = 0.107, alpha = 17.02) {
  args <- as.list(environment())

  if (!all(sapply(args, is.numeric)))
    stop("All arguments must be numeric.")

  l <- lengths(args)
  if (l[1] != 3)
    stop("d must be of length 3.")

  if (any(l[2:length(l)] > 1))
    stop("angle, ant.dens, LT, L0, wA, ln, wn, and alpha should be scalar values.")

  theta <- pi / 180 * angle
  wT <- 4.799 * angle ^ -0.5014
  LAi <- L0 - 0.5 * wA / tan(0.5 * theta)
  Dmax <- LAi * cos(0.5 * theta)
  N <- ant.dens * (LT + 4 * LAi)

  b <- 2 * tan(0.5 * theta) * d
  nb <- wT / (ln * wn) * (1 - wT * tan(0.5 * theta)) * b ^ 2
  f <- LT + (1 - d[1] / Dmax) * LAi + (1 - d[3] / Dmax) * LAi +
    mid_dist(d[1], d[2], theta, Dmax) + mid_dist(d[2], d[3], theta, Dmax) +
    b[1] + b[2] + b[3]

  (N - nb[1] / alpha - nb[2] / alpha - nb[3] / alpha) / f
}
