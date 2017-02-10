mid_dist <- function(d1, d2, theta, Dmax){

  return(abs((1/cos(0.5*theta))*(Dmax - d1 - d2)))

}
