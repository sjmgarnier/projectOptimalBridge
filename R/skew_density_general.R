#' @title Density Function for Asymmetric Apparatus
#'
#' @description This function models the density as a function of
#'              distance for the asymmetric appratus. 
#' 
#' @param d vector of distances values 
#'       
#' @param angle1 Angle left arm makes with vertical (in degrees)
#' 
#' @param angle2 Angle right arm makes with vertical (in degrees)
#'   
#' @param ant.dens Density of ants, default value =  2.2
#' 
#' @param LT Overall trail length, default value = 100 (cm)
#' 
#' @param LS0 Length of left arm, default value = 22 (cm)
#' 
#' @param LL0 Length of right arm, default value = 44 (cm)
#' 
#' @param wA Width of apparatus, default value = 3.3 (cm)
#' 
#' @param ln Length of avg ant default value = 0.691 (cm)
#' 
#' @param wn  Width of avg ant default value = 0.107 (cm)
#' 
#' @param alpha Free fitting parameter value from Reid et al, default value 17.02
#' 
#' @return Density \code{rho}
#'
#' @author Jason Graham, \email{jason.graham@scranton.edu}
#'
#'  @export
#'    
skew_density_general <- function(d,
                                   angle1,
                                   angle2,
                                   ant.dens = 2.2,
                                   LT = 100,
                                   LS0 = 22,
                                   LL0 = 44,
                                   wA = 3.3,
                                   ln = 0.691,
                                   wn = 0.107,
                                   alpha = 17.02 ){

   # Function defining density function for skew-bridge model
  theta <- pi/180*angle1   # angles in radians
  phi <- pi/180*angle2     # angles in radians
  wT <-  4.799*(angle1 + angle2)^(-0.5014)
  LAS <- LS0 - 0.5*wA/tan(0.5*(theta + phi)) # (cm)
  LAL <- LL0 - 0.5*wA/tan(0.5*(theta + phi)) # (cm)
  LA <- LAS + LAL
  N <- ant.dens*(LT+LA)
  Dmax1 <- LAS*cos(theta) # maximum left distance from appex (cm)
  Dmax2 <- LAL*cos(phi)   # maximum right distance from appex (cm)


  ### specification of function ###

  # length of bridge squared (from law of cosines)
  b.squared <- (LAS/Dmax1*d[1])^2 + (LAL/Dmax2*d[2])^2 - 2*d[1]*d[2]*((LAS*LAL)/(Dmax1*Dmax2))*cos(theta+phi)
  b <- sqrt(b.squared)  # length of living bridge (cm)
  f <- LT + LA + b - (LAS/Dmax1)*d[1] - (LAL/Dmax2)*d[2] # length of path of traval (cm)
  nb <- (wT/(ln*wn))*(1 - wT*tan(0.5*(theta+phi)))*b.squared # number of ants sequesterd for bridge

  rho <- (N - nb/alpha)/f # density function

  return(rho)
}
