### parameters ###
#ant.dens <- 2.2
#LT <- 100 # overall trail length (cm)
#LS0 <- 22 # length of left arm (cm)
#LL0 <- 44 # length of right arm (cm)
#wA <- 3.3 # width of apparatus (cm)
#ln <- 0.691 # length of avg ant (cm)
#wn <- 0.107 # width of avg ant (cm)
#wT <- 4.799*(angle1 + angle2)^(-0.5014)  # ratio of bridge width/length (value from Reid et al)
#alpha <- 17.02 #free fitting parameter value from Reid et al


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
