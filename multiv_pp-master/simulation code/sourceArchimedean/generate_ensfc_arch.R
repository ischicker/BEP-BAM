# generate multivariate ensemble forecasts 
# Note that parts of the notation may differ from the notation in the paper
#   ... in particular regarding the numbering of the settings which was changed during the revision (2 -> S1; 3 -> 2; 4 -> 3A, new: 3B)

# input:
#   model: character string; indicating which model is used for generating the forecasts
#   nout: number of multivariate observations to be generated as evaluation period
#   ninit: additional initial training period for model estimation purposes
#   nmembers: number of ensemble members
#   d: dimension of the multivariate vectors
#   ... additional parameters, depending on the chosen model

# output:
#   list of two arrays "ensfc_init" and "ensfc" containing the ensemble forecasts
#   dimensions: ninit, nmembers, d; and nout, nmembers, d
#   content in array dimensions:
#     first dimension: forecast instance 
#     second dimension: ensemble member
#     third dimension: dimension in multivariate setting

# Example calls:
# generate_ensfc(1, 10, 10, 10, 3, theta=5, copula= "Frank")


generate_ensfc <- function(model, nout, ninit, nmembers, ...){
  require(MASS)
  
  d <- list(...)$d
  
  # check input 
  if(any(!is.numeric(c(nout, ninit, nmembers, d)))){
    stop("Input 'nout', 'ninit', 'nmembers' and 'd' need to be numeric of length 1")
  }
  m <- nmembers
  
  # initialize output arrays
  ensfc_init <- array(NA, dim = c(ninit, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  
  # Setting 1
  if(model == 1){
    
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("theta","copula")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # Assign dependence parameter and copula from input
    theta <- input$theta
    copula <- input$copula
    
    if (copula == "Frank"){
      # -infty < theta < infty and theta != 0
      cop <- frankCopula(param = theta, dim = d)
    } else if (copula == "Gumbel") {
      # 1 <= theta < infty
      cop <- gumbelCopula(param = theta, dim = d)
    } else if (copula == "Clayton") {
      # -1 <= theta < infty and theta != 0
      cop <- claytonCopula(param = theta, dim = d)
    } else {
      stop(paste("Incorrect value for 'copula' input"))
    }
    
    
    paramMargins <- list()
    
    for (i in 1:d){
      paramMargins[[i]] <- list(mean = 0, sd = 1)
    }
    
    # Generate observations with standard normal marginals
    mvDistribution <- mvdc(copula=cop, margins=rep("norm", d),
                           paramMargins=paramMargins)
    
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- rMvdc(m, mvDistribution)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- rMvdc(m, mvDistribution)
      ensfc[nn,,] <- tmp
    }
  }
  
  # Setting 2 (Sampled tau)
  if(model == 2){
    
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("tau","copula")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # Assign dependence parameter and copula from input
    tau <- input$tau
    copula <- input$copula
    
    if (copula == "Frank"){
      # -infty < theta < infty and theta != 0
      cop <- frankCopula(param = tauToParam("Frank",tau), dim = d)
    } else if (copula == "Gumbel") {
      # 1 <= theta < infty
      cop <- gumbelCopula(param = tauToParam("Gumbel",tau), dim = d)
    } else if (copula == "Clayton") {
      # -1 <= theta < infty and theta != 0
      cop <- claytonCopula(param = tauToParam("Clayton",tau), dim = d)
    } else {
      stop(paste("Incorrect value for 'copula' input"))
    }
    
    
    paramMargins <- list()
    
    for (i in 1:d){
      paramMargins[[i]] <- list(mean = 0, sd = 1)
    }
    
    # Generate observations with standard normal marginals
    mvDistribution <- mvdc(copula=cop, margins=rep("norm", d),
                           paramMargins=paramMargins)
    
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- rMvdc(m, mvDistribution)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- rMvdc(m, mvDistribution)
      ensfc[nn,,] <- tmp
    }
  }
  
  if (model == 3){
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("eps", "sigma", "rho")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # assign model-specific parameters from input
    eps <- input$eps
    sigma <- input$sigma
    rho <- input$rho
    
    # correlation matrix
    S <- matrix(NA, d, d)
    for(i in 1:d){
      for(j in 1:d){
        S[i,j] <- sigma*rho^(abs(i-j))
      }
    }
    
    # mean vector
    mu <- rep(eps, d)
    
    # generate forecasts
    for(nn in 1:ninit){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc_init[nn,,] <- tmp
    }
    for(nn in 1:nout){
      tmp <- mvrnorm(n = m, mu = mu, Sigma = S)
      ensfc[nn,,] <- tmp
    }
  }

  return(list("ensfc_init" = ensfc_init, "ensfc" = ensfc))
}
