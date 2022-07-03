# File that contains the settings for each simulation run


# Model 1 : Standard Gaussian Marginals for fixed theta
# Model 2 is for a sampled Kendall's tau from which the copula parameters are estimated
# Model 3 is a reproduction of results Lerch with also Archimedean copula fitting

getModelSettings <- function(modelSetting) {
  .GlobalEnv$observationsModel <- modelSetting
  .GlobalEnv$forecastModel <- modelSetting
  .GlobalEnv$modelSetting <- modelSetting
  
  # Parameter that influences reliability of methods such as Ssh and ECC.S --> large value for MC reps also does the trick
  .GlobalEnv$randomRepetitions <- 1
  
  # Save data
  .GlobalEnv$Rdata_dir <- "../Data/Rdata/Archimedean" # directory to save Rdata files to
  .GlobalEnv$Rout_dir <- "../Data/Rout/Archimedean" # directory to save Rout files to
  
  
  ### Mismatch in parameter value (theta and theta0) ###
  if (modelSetting == 1) {
    # Dependence parameter between weather variables for observations
    .GlobalEnv$input_theta0 <- c(5, 10)
    # Dependence parameter between weather variables for ensemble forecasts
    .GlobalEnv$input_theta <- c(5, 10)
    # Copula type
    .GlobalEnv$input_copula <- c("Frank","Gumbel", "Clayton")
    # Dimension of multi-index l (weather variable, position and look-ahead time)
    .GlobalEnv$input_d <- 3
    # MC steps for boxplots
    .GlobalEnv$MC_reps <- 100
    
    # Put the parameters in a grid
    .GlobalEnv$input_par <- expand.grid(input_theta0, input_theta, input_copula, input_d)
    names(.GlobalEnv$input_par) <- c("theta0", "theta", "copula",  "d")
    
    # Simulation parameters
    .GlobalEnv$evalDays <- 150
    .GlobalEnv$timeWindow <- 30
    .GlobalEnv$trainingDays <- 50
    .GlobalEnv$ensembleMembers <- 50
    
  } 
  ### Sample Kendall's tau repetitions times ###
  else if (modelSetting == 2) {
    # Copula type
    .GlobalEnv$input_copula <- c("Frank","Gumbel", "Clayton")
    # Dimension of multi-index l (weather variable, position and look-ahead time)
    .GlobalEnv$input_d <- 3
    
    # Repetitions
    .GlobalEnv$repetitions <- 10
    
    # Put the parameters in a grid
    .GlobalEnv$input_par <- expand.grid(input_copula, input_d, 1:repetitions)
    names(.GlobalEnv$input_par) <- c("copula",  "d", "repetition")
    
    # MC steps for boxplots
    .GlobalEnv$MC_reps <- 100
    
    # Training days should be >= m for Schaake shuffle
    .GlobalEnv$evalDays <- 150
    .GlobalEnv$timeWindow <- 30
    .GlobalEnv$trainingDays <- 50
    .GlobalEnv$ensembleMembers <- 50
  } 
  ### Reproducing results Lerch ###
  else if (modelSetting == 3) {
    # Reproduce Lerch Setting 1 but with time window


    # Dimension of multi-index l (weather variable, position and look-ahead time)
    .GlobalEnv$input_d <- 3

    # Repetitions
    .GlobalEnv$repetitions <- 1

    # Fixed sd, mu0 and eps
    .GlobalEnv$eps <- 1
    .GlobalEnv$sigma <- sqrt(5)

    # Varying rho
    .GlobalEnv$rho0 <- c(0.25, 0.5, 0.75)
    .GlobalEnv$rho <- c(0.25, 0.5, 0.75)

    # Put the parameters in a grid
    .GlobalEnv$input_par <- expand.grid(rho0, eps, sigma, rho, input_d, 1:repetitions)
    names(.GlobalEnv$input_par) <- c("rho0", "eps", "sigma", "rho", "d", "repetition")

    # MC steps for boxplots
    .GlobalEnv$MC_reps <- 100

    # Training days should be >= m for Schaake shuffle
    .GlobalEnv$evalDays <- 150
    .GlobalEnv$timeWindow <- 30
    .GlobalEnv$trainingDays <- 50
    .GlobalEnv$ensembleMembers <- 50

  }
}

