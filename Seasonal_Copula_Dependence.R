# rm(list=ls())

# library(ggplot2)
# library(ggpubr)
# library(purrr)
# library(qqplotr)
# require(MASS)
# library(copula)
# library("hydroTSM")
# 
# setwd("C:/Users/20192042/OneDrive - TU Eindhoven/University/Bachelor/Year 3/BEP - BAM/Code/multiv_pp-master/simulation code")
# 
# # "source" files for functions
# dir <- "./sourceArchimedean/"
# source(paste0(dir, "postprocess_ensfc_arch.R"))
# # source(paste0(dir, "mvpp_arch.R"))
# # source(paste0(dir, "evaluation_functions_arch.R"))
# source(paste0(dir, "CopulaParameter.R"))
# source("ECC_T2M_Emos_subfunctions.R")


# Load data
# source("getData.R")
# getData()



# columnNames <- c("Group", "Season", "Clayton", "Frank", "Gumbel","Max")
# likelihoodList <- data.frame(matrix(ncol = length(columnNames), nrow = 0))
# colnames(likelihoodList) <- columnNames
# 
# 
# 
# # Get seasons
# data1$season <- time2season(data1$time, out.fmt = "seasons")
# data2$season <- time2season(data2$time, out.fmt = "seasons")
# data3$season <- time2season(data3$time, out.fmt = "seasons")



copulaPerSeason <- function(data, season, group) {
  

  # Univariate Postprocessing - UVPP
  
  
  # Stations and days from the data
  stations <- unique(data$stat)
  days <- sort(unique(data$td))
  trainingDays <- length(days)
  
  
  # 1 weather variable over 1 lead time with d stations
  d <- length(stations)
  
  # Ensemble members
  m <- sum(grepl("laef", names(data)))
  ensembleMembers <- sapply(1:m, FUN = function(x) paste0("laef", x))
  
  # Create ensfc and obs data structure
  obs_init <- array(NA, dim = c(trainingDays, d))
  obs <- array(NA, dim = c(d))
  
  # Dimension only consists of stations
  for (station in stations) {
    # Element of 1..d
    index <- match(station,stations)
    
    # Reformat days to be in 1..days
    for (day in 1:length(days)) {
      # Extract forecast for day and dim = index
      tryCatch(
        {
          dat <- subset(data, td == days[day] & stat == station)
          if (day <= trainingDays) {
            obs_init[day, index] <- unlist(dat["obs"], use.names = FALSE)
          } else{
            obs[day - trainingDays, index] <- unlist(dat["obs"], use.names = FALSE)
          }
          
        }
        ,
        error = function(cond) {
          print(cond)
        }
      )
    }
  }
  
  uvpp <- list()
  
  for (station in stations) {
    
    statIndex <- match(station,stations)
    temp <- subset(data, stat == station)
    
    # EMOS over entire period
    emos.result <- emos_T2M_mean_singleForecast(temp, dim(temp)[1]-1)
    emos.result$stat <- statIndex
    
    uvpp <- rbind(uvpp, emos.result)
    
    
  }
  
  .GlobalEnv$obs_init <- obs_init
  .GlobalEnv$obs <- obs
  .GlobalEnv$uvpp <- uvpp


  # d <- dim(obs_init)[2]
  # # The input for the copulas are the CDF values of the margins
  # obs_train <- obs_all
  # obs_train_CDF <- c()
  # mean_vector <- c()
  # sd_vector <- c()
  # for (dd in 1:d) {
  #   # Use EMOS values for marginals
  # 
  #   dat <- subset(uvpp, stat == dd)
  #   averagedMean <- unlist(unname(dat["ens_mu"]))
  #   # print(averagedMean)
  #   averagedSd <- unlist(unname(dat["ens_sd"]))
  #   # print(averagedSd)
  # 
  #   
  #   obs_train_CDF <- cbind(obs_train_CDF, pnorm(obs_train[,dd], mean = averagedMean, sd = averagedSd))
  #   
  #   # Add for later use
  #   mean_vector <- c(mean_vector, averagedMean)
  #   sd_vector <- c(sd_vector, averagedSd)
  # }
  # 
  # 
  # try({
  #   clayton <- fitCopula(claytonCopula(dim = d), data = obs_train_CDF, method="mpl", start = 2,optim.control = list(maxit=1000), upper = 100)
  #   claytonLogLik <- clayton@loglik
  #   frank <- fitCopula(frankCopula(dim = d), data = obs_train_CDF, method="mpl", start = 5,optim.control = list(maxit=1000), upper = 100)
  #   frankLogLik <- frank@loglik
  #   gumbel <- fitCopula(gumbelCopula(dim = d), data = obs_train_CDF, method="mpl", start = 2,optim.control = list(maxit=1000), upper = 100)
  #   gumbelLogLik <- gumbel@loglik
  #   
  #   copulas <- c("Clayton", "Frank", "Gumbel")
  #   index <- which.max(c(claytonLogLik, frankLogLik, gumbelLogLik))
  #   
  #   .GlobalEnv$likelihoodList[nrow(.GlobalEnv$likelihoodList) + 1,] = c(group, season, claytonLogLik, frankLogLik, gumbelLogLik, copulas[index])
  # }
  #   )
  
  

}

copulaPerSeason(subset(data1, season == "summer"), "summer", 1)



# for (season in c("summer", "autumn", "winter", "spring")) {
#   copulaPerSeason(subset(data1, season == season), season, 1)
#   copulaPerSeason(subset(data2, season == season), season, 2)
#   copulaPerSeason(subset(data3, season == season), season, 3)
# }










