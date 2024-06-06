psurv_norm <- function(x, mean=0, sd=1) {
  return(1-pnorm(x,mean=mean, sd=sd))
}

qsurv_norm <- function(x, mean=0, sd=1) {
  return(1-qnorm(x,mean=mean, sd=sd))
}

dsurv_norm <- function(x, mean=0, sd=1) {
  return(-dnorm(x,mean=mean, sd=sd))
}


library(copula)
mvpp_fixed <- function(method, ensfc, ensfc_init, obs, obs_init, postproc_out, 
                 EMOS_sample = NULL, ECC_out = NULL, timeWindow, ecc_m, 
                 uvpp = NULL, training_df){
  
  # include some checks for inputs
  
  # generate array for ouput
  params <- array(NA, dim = dim(ensfc)[1])
  chosenCopula <- array(NA, dim = dim(ensfc)[1])
  n <- dim(ensfc)[1]
  if (method %in% c("ECC", "dECC")) {
    m <- dim(ensfc)[2]
  } else {
    m <- ecc_m
  }
  d <- dim(ensfc)[3]
  mvppout <- array(NA, dim = c(n,m,d))
  
  
  
  # Schaake shuffle code
  if(method == "SSh-H"){
    
    # concatenate obs_init and obs arrays to sample from available forecast cases later on
    obs_all <- rbind(obs_init, obs)
    
    # reorder post-processed forecast sample according to past observations
    for(nn in 1:n){
      # choose set of past forecast cases to determine dependence template
      #   ... this way, a new set of IDs is drawn for every forecast instance
      #   ... this needs to depend on nn in a more suitable manner if there is temporal change in the simulation setup
      obs_IDs <- training_df$TrainingDays[[nn]]
      for(dd in 1:d){
        obs_tmp <- obs_all[obs_IDs, dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(obs_tmp, ties.method = "random")]
      }
    }
    # end of SSh code   
  }
  
  if(method == "SSh-I14"){
    
    # concatenate obs_init and obs arrays to sample from available forecast cases later on
    obs_all <- rbind(obs_init, obs)
    
    # reorder post-processed forecast sample according to past observations
    for(nn in 1:n){
      # choose set of past forecast cases to determine dependence template
      #   ... this way, a new set of IDs is drawn for every forecast instance
      #   ... this needs to depend on nn in a more suitable manner if there is temporal change in the simulation setup
      obs_IDs <- training_df$TrainingDays[[nn]]
      for(dd in 1:d){
        obs_tmp <- obs_all[obs_IDs, dd]
        mvppout[nn, , dd] <- EMOS_sample[nn, , dd][rank(obs_tmp, ties.method = "random")]
      }
    }
    # end of SSh code   
  }
  
  # GCA code
  if(method == "GCA"){
    require(MASS)
    
    # concatenate obs_init and obs arrays to determine covariance matrix for Gaussian copulas
    obs_all <- rbind(obs_init, obs)
    
    for(nn in 1:n){
      
      # # Only last measurements
      obs_IDs <- training_df$TrainingDays[[nn]]
      obs_train <- obs_all[obs_IDs, ]
      
      # Latent Gaussian observations
      obs_latent_gaussian <- c()
      mean_values <- c()
      sd_values <- c()
      for(dd in 1:d){
        
        dat <- subset(uvpp, stat == dd)
        averagedMean <- unlist(unname(dat[nn,]["ens_mu"]))
        averagedSd <- unlist(unname(dat[nn,]["ens_sd"]))
        
        # NEW
        delta = 1e-3
        obs_train_CDF_raw <- pnorm(obs_train[,dd], mean = averagedMean, sd = averagedSd)
        obs_latent_gaussian <- cbind(obs_latent_gaussian, qnorm(pmin(pmax(obs_train_CDF_raw, delta), 1 - delta)))
        
        mean_values <- c(mean_values, averagedMean)
        sd_values <- c(sd_values, averagedSd)
      }
      
      # estimate covariance matrix
      cov_obs <- cov(obs_latent_gaussian)
      # draw random sample from multivariate normal distribution with this covariance matrix
      # Make sure to get numeric values
      
      # Dependence structure by Copula
      mvsample <- mvrnorm(n = m, mu = rep(0, d), Sigma = cov_obs)
      
      # impose dependence structure on post-processed forecasts
      for(dd in 1:d){
        
        mvppout[nn, , dd] <- mvsample[, dd] * sd_values[dd] + mean_values[dd]
        
      }
    }
    .GlobalEnv$t <- mvppout
    
    # end of GCA code 
  }
  
  # Archimedean copula code
  if(any(match(c("Clayton", "Frank", "Gumbel", "Surv_Gumbel", "GOF"), method, nomatch = 0) != 0 ) ){
    require(MASS)
    
    # concatenate obs_init and obs arrays to determine covariance matrix for Gaussian copulas
    obs_all <- rbind(obs_init, obs)
    
    for(nn in 1:n){
      # Use a shifting time window of length obs_init to compute copula parameters
      obs_IDs <- training_df$TrainingDays[[nn]]
      obs_train <- obs_all[obs_IDs, ]
      
      d <- dim(obs_train)[2]
      # The input for the copulas are the CDF values of the margins
      obs_train_CDF <- c()
      mean_vector <- c()
      sd_vector <- c()
      for (dd in 1:d) {
        # Use EMOS values for marginals
        if (!is.null(uvpp)) {
          dat <- subset(uvpp, stat == dd)
          averagedMean <- unlist(unname(dat[nn,]["ens_mu"]))
          # print(averagedMean)
          averagedSd <- unlist(unname(dat[nn,]["ens_sd"]))
          # print(averagedSd)
        } else {
          par <- postproc_out[nn, dd, ]
          averagedMean <- par[1]
          averagedSd <- par[2]
        }
        
        obs_train_CDF <- cbind(obs_train_CDF, pnorm(obs_train[,dd], mean = averagedMean, sd = averagedSd))
        
        # Add for later use
        mean_vector <- c(mean_vector, averagedMean)
        sd_vector <- c(sd_vector, averagedSd)
      }
      
      # Estimate the parameter and copula - itau method does not converge for some cases (without giving a warning/ error and runs indefinitely)
      # Copula parameters are bounded to prevent generating inf samples
      if (method == "Clayton") {
        fitcop <-tryCatch({
          fitCopula(claytonCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        }, error = function(e) {
          # print(e)
          indepCopula(dim=d)
        })
      }
      else if (method == "Frank") {
        fitcop <-tryCatch({
          fitCopula(frankCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        }, error = function(e) {
          # print(e)
          indepCopula(dim=d)
        })
      }
      else if (method == "Gumbel") {
        fitcop <-tryCatch({
          fitCopula(gumbelCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        }, error = function(e) {
          # print(e)
          indepCopula(dim=d)
        })
      } else if (method == "Surv_Gumbel") {
        fitcop <-tryCatch({
          fitCopula(gumbelCopula(dim = d), data = 1-obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        }, error = function(e) {
          # print(e)
          indepCopula(dim=d)
        })
      }else {
        stop("Incorrect copula")
      }
      
      if (!(class(fitcop) == "indepCopula")){
        cop <- fitcop@copula
        chosenCopula[nn] <- method
        
        # Save the fitted parameter
        params[nn] <- fitcop@estimate
      } else {
        chosenCopula[nn] <- "indep"
        cop <- fitcop
      }
      
      # draw random sample from multivariate normal distribution with this estimated copula
      paramMargins <- list()
      surv_paramMargins <- list()
      
      for (i in 1:d){
        paramMargins[[i]] <- list(mean = mean_vector[i], sd = sd_vector[i])
        surv_paramMargins[[i]] <- list(mean = -mean_vector[i], sd = sd_vector[i])
      }
      
      # Generate observations with normal marginals
      if (method == "Surv_Gumbel") {
        mvDistribution <- mvdc(copula=cop, margins=rep("surv_norm", d),
                               paramMargins=surv_paramMargins)
      } else {
        mvDistribution <- mvdc(copula=cop, margins=rep("norm", d),
                               paramMargins=paramMargins)
      }
      
      
      # Make sure to get numeric values
      max_reps <- 2
      rep <- 0
      repeat {
        rep <- rep + 1
        # Dependence structure by Copula
        mvsample <- rMvdc(m, mvDistribution)
        
        if (all(is.finite(mvsample)) | rep >= max_reps) {
          break
        }
        
        set.seed(sample(1:1000,1))
        
        # print(rep)
        # print("Repeating...2")
      }
      
      for(dd in 1:d){
        
        # Translate CDF values back to forecasts
        mvppout[nn, , dd] <- mvsample[,dd]
      }
      
    }
    
    # end of Archimedean Copula code
    
    
  }
  
  # GCA from copula
  if(method == "CopGCA"){
    require(MASS)
    
    # concatenate obs_init and obs arrays to determine covariance matrix for Gaussian copulas
    obs_all <- rbind(obs_init, obs)
    
    for(nn in 1:n){
      # Use a shifting time window of length obs_init to compute copula parameters
      obs_train <- obs_all[(dim(obs_init)[1]+nn - timeWindow):(dim(obs_init)[1]+nn-1), ]
      
      d <- dim(obs_train)[2]
      # The input for the copulas are the CDF values of the margins
      obs_train_CDF <- c()
      mean_vector <- c()
      sd_vector <- c()
      for (dd in 1:d) {
        # Use EMOS values for marginals
        dat <- subset(uvpp, stat == dd)
        averagedMean <- unlist(unname(dat[nn,]["ens_mu"]))
        # print(averagedMean)
        averagedSd <- unlist(unname(dat[nn,]["ens_sd"]))
        # print(averagedSd)
        
        delta = 1e-3
        obs_train_CDF_raw <- pnorm(obs_train[,dd], mean = averagedMean, sd = averagedSd)
        obs_train_CDF <- cbind(obs_train_CDF, pmin(pmax(obs_train_CDF_raw, delta), 1 - delta))
        
        # Add for later use
        mean_vector <- c(mean_vector, averagedMean)
        sd_vector <- c(sd_vector, averagedSd)
      }
      
      fitcop <-tryCatch({
        fitCopula(normalCopula(dim = d, dispstr = "un"), data = obs_train_CDF, method="mpl", optim.control = list(maxit=1000))
      }, error = function(e) {
        print(e)
        indepCopula(dim=d)
      })
      
      if (!(class(fitcop) == "indepCopula")){
        cop <- fitcop@copula
      } else {
        cop <- fitcop
      }
      
      
      # draw random sample from multivariate normal distribution with this estimated copula
      paramMargins <- list()
      surv_paramMargins <- list()
      
      for (i in 1:d){
        paramMargins[[i]] <- list(mean = mean_vector[i], sd = sd_vector[i])
      }
      
      
      mvDistribution <- mvdc(copula=cop, margins=rep("norm", d),
                             paramMargins=paramMargins)
      
      
      
      # Make sure to get numeric values
      max_reps <- 2
      rep <- 0
      repeat {
        rep <- rep + 1
        # Dependence structure by Copula
        mvsample <- rMvdc(m, mvDistribution)
        
        if (all(is.finite(mvsample)) | rep >= max_reps) {
          break
        }
        
        set.seed(sample(1:1000,1))
        
        # print(rep)
        # print("Repeating...2")
      }
      
      for(dd in 1:d){
        
        # Translate CDF values back to forecasts
        mvppout[nn, , dd] <- mvsample[,dd]
      }
      
    }
    
    # end of copula gca
    
    
  }
  
  return(list("mvppout" = mvppout, "params" = params, "chosenCopula" = chosenCopula))
  # in random methods: distinguish cases with and without given EMOS_sample, maybe only handle that with sample at first, rest can be included later on
  # if no sample is given, a new one has to be generated, as done for the EMOS methods themselves
}

