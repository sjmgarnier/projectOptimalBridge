#' @title Density Function for Asymmetric Apparatus
#'
#' @description This function computes the density as a function of distance for
#'  the asymmetric appratus.
#'
#' @param d A length 2 numeric vector of distances values in cm.
#'
#' @param angle1 A scalar for the angle in deg that the left arm makes with
#'  vertical.
#'
#' @param angle2 A scalar for the angle in deg that the right arm makes with
#'  vertical.
#'
#' @param ant.dens A scalar for the density of ants in ants/cm (default: 2.2).
#'
#' @param LT A scalar for the overall trail length in cm (default: 100).
#'
#' @param LS0 A scalar for the length of the left arm in cm (default: 22).
#'
#' @param LL0 A scalar for the length of the right arm in cm (default: 44).
#'
#' @param wA A scalar for the width of the arms of the apparatus in cm (default: 3.3).
#'
#' @param ln A scalar for the average length of an ant in cm (default:0.691).
#'
#' @param wn  A scalar for the average width of an ant in cm (default: 0.107).
#'
#' @param alpha A scalar for the free fitting parameter value from Reid et al.
#'  (default: 17.02).
#'
#' @return A vector of the of the same length as \code{d}.
#'
#' @author Jason Graham, \email{jason.graham@@scranton.edu}
#'
#' @author Simon Garnier, \email{garnier@@njit.edu}
#'
#' @export
#'
skew_density_general <- function(d, angle1, angle2, ant.dens = 2.2, LT = 100,
                                 LS0 = 22, LL0 = 44, wA = 3.3, ln = 0.691,
                                 wn = 0.107, alpha = 17.02) {
  args <- as.list(environment())

  if (!all(sapply(args, is.numeric)))
    stop("All arguments must be numeric.")

  l <- lengths(args)
  if (l[1] != 2)
    stop("d must be of length 2.")

  if (any(l[2:length(l)] > 1))
    stop("angle1, angle2, ant.dens, LT, L0, wA, ln, wn, and alpha should be scalar values.")

  theta <- pi / 180 * angle1
  phi <- pi / 180 * angle2
  wT <- 4.799 * (angle1 + angle2) ^ -0.5014
  LAS <- LS0 - 0.5 * wA / tan(0.5 * (theta + phi))
  LAL <- LL0 - 0.5 * wA / tan(0.5 * (theta + phi))
  LA <- LAS + LAL
  N <- ant.dens * (LT + LA)
  Dmax1 <- LAS * cos(theta)
  Dmax2 <- LAL * cos(phi)

  b.squared <- (LAS / Dmax1 * d[1]) ^ 2 + (LAL / Dmax2 * d[2]) ^ 2 -
    2 * d[1] * d[2] * ((LAS * LAL) / (Dmax1 * Dmax2)) * cos(theta + phi)
  b <- sqrt(b.squared)
  f <- LT + LA + b - (LAS / Dmax1) * d[1] - (LAL / Dmax2) * d[2]
  nb <- (wT / (ln * wn)) * (1 - wT * tan(0.5 * (theta + phi))) * b.squared

  (N - nb / alpha) / f
}
