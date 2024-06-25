#####################################################################
#  Similarity criterion as in Schefzik 2016 (Sim Schaake) Formula (4)
#  Elisa Perrone
#  2021 - 05 - 25
####################################################################

Sim <- function(t,td, uvpp, ensfc){
  
  # Ensemble members and number of stations
  m <- dim(ensfc)[2]
  d <- dim(ensfc)[3]
  
  Mean<-0
  Sd<-0
  
  for(dd in 1:d){
    
    # Retrieve the EMOS statistics for this station
    dat <- subset(uvpp, stat == dd)
    
    # Extract the variables
    xbar_t  <- unlist(unname(dat[t, ]["ens_mu"]))
    s_t     <- unlist(unname(dat[t, ]["ens_sd"]))
    
    xbar_td <- unlist(unname(dat[td, ]["ens_mu"]))
    s_td    <- unlist(unname(dat[td, ]["ens_sd"]))
    
    # Use the statistics to compute the criterion
    Mean <- Mean + (xbar_t - xbar_td)^2
    Sd   <- Sd   + (s_t - s_td)^2
  }
  
  z <- sqrt( (1/d) * sum(Mean) + (1/d) * sum(Sd)  )
  return(z)
}