# code to generate multivariate observations 
# Note that parts of the notation may differ from the notation in the paper
#   ... in particular regarding the numbering of the settings which was changed during the revision (2 -> S1; 3 -> 2; 4 -> 3A, new: 3B)

# input:
#   model: numeric, indicating which model is used for the observations
#   nout: number of multivariate observations to be generated as evaluation period 
#   ninit: additional initial training period for model estimation purposes
#   d: dimension of the multivariate vectors
#   ... additional parameters, depending on the chosen model

# output:
#   list of two arrays "obs_init" and "obs" containing the observations
#   dimensions: ninit, d; and nout; d
#   first dimension: forecast instance; second dimension: dimension in multivariate setting


# Example calls:
# generate_obs(1, 10, 10, 3, theta0=0.5, copula= "Frank")
# example_plot(5, "Clayton")
# example_plot(5, "Frank")
# example_plot(5, "Gumbel")

generate_obs <- function(model, nout, ninit, d, ...){
  require(MASS)
  require(copula)
  require(scatterplot3d)
  
  # check input 
  if(any(!is.numeric(c(nout, ninit, d)))){
    stop("Input 'nout', 'ninit' and 'd' need to be numeric of length 1")
  }
  
  # initialize output arrays
  # Obs_init denotes the observations for the training period and obs for the evaluation period
  obs_init <- array(NA, dim = c(ninit, d))
  obs <- array(NA, dim = c(nout, d))
  
  
  # Setting 1 (Gaussian marginals)
  if(model == 1){
    
    # check if appropriate additional parameters are given
    input <- list(...)
    required <- c("theta0","copula")
    ind_check <- match(required, names(input), nomatch = 0)
    if(any(ind_check == 0)){
      stop(paste("Missing additional model-specific parameter",
                 paste("Given input:", paste(names(input), collapse=", ")),
                 paste("Required input:", paste(required, collapse=", ")),
                 sep="\n")
      )
    }
    
    # Assign dependence parameter and copula from input
    theta0 <- input$theta0
    copula <- input$copula
    
    if (copula == "Frank"){
      # -infty < theta0 < infty and theta0 != 0
      cop <- frankCopula(param = theta0, dim = d)
    } else if (copula == "Gumbel") {
      # 1 <= theta0 < infty
      cop <- gumbelCopula(param = theta0, dim = d)
    } else if (copula == "Clayton") {
      # -1 <= theta0 < infty and theta0 != 0
      cop <- claytonCopula(param = theta0, dim = d)
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
    
    # Random sample from the mvdc
    obs_init <- rMvdc(ninit, mvDistribution)
    obs <- rMvdc(nout, mvDistribution)
    
    
  }
  

  
  return(list("obs_init" = obs_init, "obs" = obs))
}

?gumbelCopula

example_plot <- function(theta0, copula) {
  library("car")
  library("rgl")
  data <- generate_obs(1, 10000, 1, 3, theta0=theta0, copula= copula)$obs
  par3d(windowRect = c(20, 30, 800, 800))
  scatter3d(data[,1],data[,2], data[,3],point.col = "blue", surface=FALSE, size=40)
}
