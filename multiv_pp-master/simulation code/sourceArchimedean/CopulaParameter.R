# Mapping between Kendall's tau and parameter

require(FactorCopula)

paramToTau <- function(copula, param) {
  if (copula == "Clayton") {
    return(param/(param + 2))
  } else if (copula == "Frank") {
    return(par2tau("frk", param))
  } else if (copula == "Gumbel") {
    return(par2tau("gum", param))
  } else {
    stop("Not a suitable copula")
  }
}

tauToParam <- function(copula, tau) {
  if (copula == "Clayton") {
    return(2 * tau / (1-tau))
  } else if (copula == "Frank") {
    return(tau2par("frk", tau))
  } else if (copula == "Gumbel") {
    return(tau2par("gum", tau))
  } else {
    stop("Not a suitable copula")
  }
}