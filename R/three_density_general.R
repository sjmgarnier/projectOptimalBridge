### parameters ###
#ant.dens <- 2.2
#LT <- 100 # overall trail length (cm)
#ln <- 0.691 # length of avg ant (cm)
#wn <- 0.107 # width of avg ant (cm)
#wA <- 3.3 # width of appratus (cm)
#alpha <- 17.02 #free fitting parameter from Reid et al


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
