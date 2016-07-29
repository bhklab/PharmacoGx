.GR <- function(x, pars, tau) {
  #GR takes in a vector of log concentrations, a vector of DRC parameters from the .Hill()
  #function, and a coefficient tau equal to the number of doubling times occuring between
  #the start of the experiment and the taking of the viability measurements. It then returns
  #the GR-value associated with those conditions.
  return((.Hill(x, pars)) ^ (1 / tau))
}