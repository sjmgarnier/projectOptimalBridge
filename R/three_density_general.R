#' @title Density Function for Three-apex Apparatus
#'
#' @description This function models the density as a function of
#'              distance for the two-apex appratus. 
#' 
#' @param d vector of distances values 
#'       
#' @param angle Angle of apex (in degrees)
#' 
#' @param ant.dens Density of ants, default value =  2.2
#' 
#' @param LT Overall trail length, default value = 100 (cm)
#' 
#' @param L0 Length of one arm, default value = 22 (cm)
#'
#' @param ln Length of avg ant default value = 0.691 (cm)
#' 
#' @param wn  Width of avg ant default value = 0.107 (cm)
#' 
#' @param wA Width of apparatus, default value = 3.3 (cm)
#'
#' @param alpha Free fitting parameter value from Reid et al, default value 17.02
#' 
#' @return Density \code{rho}
#'
#' @author Jason Graham, \email{jason.graham@scranton.edu}
#'
#'  @export
#'    
three_density_general <- function(d,
                                  angle,
                                  ant.dens = 2.2,
                                  LT = 100,
                                  L0 = 22,
                                  ln = 0.691,
                                  wn = 0.107,
                                  wA = 3.3,
                                  alpha = 17.02 ){
  # Function specifying density for three-bridge model


  ### parameters ###
  theta <- pi/180*angle  # angle in radians
  wT <- 4.799*(angle)^(-0.5014) # ratio of bridge width/length (value from Reid et al)
  LAi <- L0 - 0.5*wA/tan(0.5*theta) # (cm)
  Dmax <- LAi*cos(0.5*theta)
  N <- ant.dens*(LT+4*LAi)

  ### specification of function ###
  ### specification of function ###
  b <- 2*tan(0.5*theta)*d
  nb <- wT/(ln*wn)*(1 - wT*tan(0.5*theta))*b^2
  f <- LT + (1 - d[1]/Dmax)*LAi + (1 - d[3]/Dmax)*LAi + mid_dist(d[1],d[2],theta,Dmax) + mid_dist(d[2],d[3],theta,Dmax) + b[1] + b[2] + b[3]

  rho <- (N - nb[1]/alpha - nb[2]/alpha - nb[3]/alpha)/f  # density function

  return(rho)
}
