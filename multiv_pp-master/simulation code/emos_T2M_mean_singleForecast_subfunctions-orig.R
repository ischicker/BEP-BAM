# -------------------------------------------------------------------
# - NAME:        emos_T2M_mean_singleForecast_subfunctions.R 
# - AUTHOR:      Moritz Lang
# - DATE:        2016-07-20
# -------------------------------------------------------------------
# - DESCRIPTION: contains subfunctions for emos_T2M_mean_singleForecast.R 
# -------------------------------------------------------------------
# - CONTENT:  - crps.normal()
#             - crps.ensemble()
#             - crps.norm()
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# CRPS calculation as in R-package 
# (adapted from Gneiting et al. 2005)
# -------------------------------------------------------------------
crps.normal <- function(param, t.obs, t.ens, t.var, K) {
   param[-1] <- param[-1]^2
   MEAN <- t.ens %*% param[1:K]
   VAR <- param[K + 1] + param[K + 2] * t.var
   SIGMA <- sqrt(VAR)
   obs.stdz <- (t.obs - MEAN) / SIGMA
   crps.obs <- SIGMA * (obs.stdz * (2 * .Internal(pnorm(obs.stdz, 0, 1, TRUE, FALSE)) - 1) + 2 * .Internal(dnorm(obs.stdz, 0, 1, FALSE)) - 1 / sqrt(pi))

   return(.Primitive("sum")(crps.obs)) # mean to sum changed!!
}
library("compiler")
crps.normal <- cmpfun(crps.normal)

# -------------------------------------------------------------------
# Continuous Rank Probability Score: on discret number of forecasts, ensemble
# (adapted from Reto Stauffer 2016)
# -------------------------------------------------------------------
crps.ensemble <- function(obs, fcst, plot = FALSE) {
   # Sorting ensemble
   fcst <- sort(fcst)
   # Ensemble empirical distribution
   fun.ecdf <- ecdf(fcst)
   # Differences and corresponding cdf's
   diff <- diff(sort(c(fcst, obs)))
   cdfs <- fun.ecdf(sort(c(fcst, obs)))[1:length(fcst)]
   # Compute crps
   crps <- sum((cdfs - (fcst > obs))^2 * diff)
   # Demo plot
   if (plot) {
      plot(fcst, fun.ecdf(fcst), type = "s", main = sprintf(" RPS score: %.3f", crps))
      abline(v = obs, col = 2)
   }
   crps
}

# -------------------------------------------------------------------
# CRPS for normal normal distribution
# (adapted from Gneiting et al. 2005)
# -------------------------------------------------------------------
crps.norm <- function(obs, mu, sd) {
   sd * ((obs - mu) / sd * (2 * pnorm((obs - mu) / sd) - 1) + 2 * dnorm((obs - mu) / sd) - 1 / sqrt(pi))
}
