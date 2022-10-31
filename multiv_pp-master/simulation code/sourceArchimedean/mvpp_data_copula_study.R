require(MASS)
library("hydroTSM")
library(copula)

seasonal_copula <- function(data, ensfc, obs, obs_init, postproc_out, ecc_m) {
  
  # Put the seasons in the dataframe
  data$season <- time2season(data$time, out.fmt = "seasons")
  
  # generate array for output
  n <- dim(ensfc)[1]
  # m <- dim(mvppout)[2]
  m <- ecc_m
  d <- dim(ensfc)[3]
  mvppout <- array(NA, dim = c(n,m,d))
  days <- sort(unique(data$td))
  
  # concatenate obs_init and obs arrays to determine covariance matrix for Gaussian copulas
  obs_all <- rbind(obs_init, obs)
  data1$season <- time2season(data1$time, out.fmt = "seasons")
  
  # Construct list of season values (note that this is the same over stations)
  temp_stat <- data[1,]$stat
  
  seasons <- c()

  # Reformat days to be in 1..days
  for (day in 1:length(days)) {
      dat <- subset(data, td == days[day] & stat == temp_stat)
      seasons <- c(seasons, time2season(dat$time, out.fmt = "seasons"))
  }
  
  columnNames <- c("Year", "Season", "Clayton_LLH","Clayton_P", "Frank_LLH", "Frank_P", "Gumbel_LLH","Gumbel_P", "Max","Max_P")
  copula_df <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
  colnames(copula_df) <- columnNames
     
  #############################   
  ## Initial training period ##
  #############################
  
  seasons_init <- seasons[1:dim(obs_init)[1]]

  for (season in c("summer", "autumm", "winter", "spring")) {
    clayton <- fitCopula(claytonCopula(dim = d), data = pobs(obs_init[seasons_init==season,]), method="mpl", start = 2,optim.control = list(maxit=1000), upper = 100)
    claytonLogLik <- clayton@loglik
    frank <- fitCopula(frankCopula(dim = d), data = pobs(obs_init[seasons_init==season,]), method="mpl", start = 5,optim.control = list(maxit=1000), upper = 100)
    frankLogLik <- frank@loglik
    gumbel <- fitCopula(gumbelCopula(dim = d), data = pobs(obs_init[seasons_init==season,]), method="mpl", start = 2,optim.control = list(maxit=1000), upper = 100)
    gumbelLogLik <- gumbel@loglik

    copulas <- c("Clayton", "Frank", "Gumbel")
    parameters <- c(clayton@copula@parameters, frank@copula@parameters, gumbel@copula@parameters)
    index <- which.max(c(claytonLogLik, frankLogLik, gumbelLogLik))

    copula_df[nrow(copula_df) + 1,] = c(0, season, 
                                       claytonLogLik, clayton@copula@parameters, 
                                       frankLogLik, frank@copula@parameters,
                                       gumbelLogLik, gumbel@copula@parameters,
                                       copulas[index], parameters[index])
    
  }
  
  #############################   
  ##    Evaluation period    ##
  #############################
  
  # Split evaluation period into discrete periods of one year
  n <- dim(obs)[1]
  year_list <- (1:n)%/%365
  
  for (year in 0:(n%/%365)) {
    # print(paste0("Year: ",year))
  
    #############################   
    ##     Use copula fits     ##
    #############################
    for (nn in (365*year+1):min(n,365*(year+1))) {
      # print(paste0("Day: ", nn))
      
      # Determine time period of current day
      season <- seasons[dim(obs_init)[1]+nn]
    
      # Extract UVPP parameters for the current day
      mean_vector <- c()
      sd_vector <- c()
      for (dd in 1:d) {
          par <- postproc_out[nn, dd, ]
          averagedMean <- par[1]
          averagedSd <- par[2]
          # Add for later use
          mean_vector <- c(mean_vector, averagedMean)
          sd_vector <- c(sd_vector, averagedSd)
      }
      
      
      # Put the parameters in the correct format
      paramMargins <- list()
      
      for (i in 1:d){
        paramMargins[[i]] <- list(mean = mean_vector[i], sd = sd_vector[i])
      }
      
  
      
      # Select proper copula
      copula_info <- copula_df[copula_df$Year == year & copula_df$Season == season,]
      copula_type <- copula_info$Max
      copula_par <- as.numeric(copula_info$Max_P)
      
      if (copula_type == "Clayton") {
        cop <- claytonCopula(dim = d, param = copula_par)
      } else if (copula_type == "Frank") {
        cop <- frankCopula(dim = d, param = copula_par)
      } else if (copula_type == "Gumbel") {
        cop <- gumbelCopula(dim = d, param = copula_par)
      }
      
      # Create the multivariate distribution and sample
      mvDistribution <- mvdc(copula=cop, margins=rep("norm", d),
                               paramMargins=paramMargins)
      mvsample <- rMvdc(m, mvDistribution)
        
      # Store the result
      for(dd in 1:d){
        
        # Translate CDF values back to forecasts
        mvppout[nn, , dd] <- mvsample[,dd]
      }
    } # END day loop
    
    
    ############################   
    ##   Update copula fits   ##
    ############################
    # Secondly: Update copula parameters based on observations of current year
    
    
    for (season in c("summer", "autumm", "winter", "spring")) {
      
      obs_past <- obs[1:min(((year+1)*365),n),]
      obs_all <- rbind(obs_init, obs_past)
      
      seasons_clipped <- seasons[1:dim(obs_all)[1]]
      
      clayton <- fitCopula(claytonCopula(dim = d), data = pobs(obs_all[seasons_clipped==season,]), method="mpl", start = 2,optim.control = list(maxit=1000), upper = 100)
      claytonLogLik <- clayton@loglik
      frank <- fitCopula(frankCopula(dim = d), data = pobs(obs_all[seasons_clipped==season,]), method="mpl", start = 5,optim.control = list(maxit=1000), upper = 100)
      frankLogLik <- frank@loglik
      gumbel <- fitCopula(gumbelCopula(dim = d), data = pobs(obs_all[seasons_clipped==season,]), method="mpl", start = 2,optim.control = list(maxit=1000), upper = 100)
      gumbelLogLik <- gumbel@loglik
      
      
      
      copulas <- c("Clayton", "Frank", "Gumbel")
      parameters <- c(clayton@copula@parameters, frank@copula@parameters, gumbel@copula@parameters)
      index <- which.max(c(claytonLogLik, frankLogLik, gumbelLogLik))
      
      
      copula_df[nrow(copula_df) + 1,] = c(year+1, season, 
                                         claytonLogLik, clayton@copula@parameters, 
                                         frankLogLik, frank@copula@parameters,
                                         gumbelLogLik, gumbel@copula@parameters,
                                         copulas[index], parameters[index])
    }
  } # END year loop
  
  return(list("mvppout" = mvppout, "copula_df" = copula_df))
}