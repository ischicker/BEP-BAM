rm(list=ls())

library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/simulation code")

# "source" files for functions 
dir <- "./sourceArchimedean/"
source(paste0(dir, "postprocess_ensfc_arch.R"))
source(paste0(dir, "mvpp_arch.R"))
source(paste0(dir, "evaluation_functions_arch.R"))
source(paste0(dir, "CopulaParameter.R"))
source("ECC_T2M_Emos_subfunctions.R")

# Load data
fName <- "INPUT-DATA_temp_2013071000-2016033000"
load(paste0("../../Data/", fName, ".Rdata")) # loads data in "input" variable

### Input consists of
# 72 different lead times
# 216 stations (input$leadtimen consists of the stations)
# 17 ensemble members (laef1, ..., leaf17)
# input$leadtimen$lt contains the lead time n
# 990 different time stamps (unique values of "initial")
# td contains the discretized time values
# initial is the time of measurement
# ytime = initial + 1
# inca (Integrated Nowcasting through Comprehensive Analyses) contains the observations
# cosmo (Consortium for Small-scale Modeling) another forecast

# Focus on one specific lead time
leadTime <- 12
data <- input[[paste0("leadtime", leadTime)]]

# Time format is YYYYMMDD00 and is converted to Date object
data$time <- as.Date(sapply(data$initial, FUN = function(x) substr(toString(x), start = 1, stop = 8)), "%Y%m%d")

#vector of all distinct stations
statNR <-  unique(input$leadtime1$stat) 


# Groups in Perrone et al. 2020 are:
# group 1: stations 188 - 189 - 191
# group 2: stations 64 - 169 - 188
# group 3: stations 19 - 150 - 172
group1 <- c(188, 189, 191)
group2 <- c(64, 169, 188)
group3 <- c(19, 150, 172)

# Select the data for those groups
data1_comp <- subset(data, stat %in% statNR[group1])
data2_comp <- subset(data, stat %in% statNR[group2])
data3_comp <- subset(data, stat %in% statNR[group3])

# Remove unnecessary data from memory
rm(input, data)

# Remove first part of data
data1 <- subset(data1_comp, initial >= 2014010100)
data2 <- subset(data2_comp, initial >= 2014010100)
data3 <- subset(data3_comp, initial >= 2014010100)

# Remove entries with NA values and rename one column
data1$obs <- data1$inca
data2$obs <- data2$inca
data3$obs <- data3$inca

data1 <- data1[!unname(apply(data1, MARGIN = 1, FUN = function(x) any(is.na(x)))),]
data2 <- data1[!unname(apply(data1, MARGIN = 1, FUN = function(x) any(is.na(x)))),]
data3 <- data1[!unname(apply(data1, MARGIN = 1, FUN = function(x) any(is.na(x)))),]


# Setting parameter for different runs with same parameters
# setting <- 1


# eval_all_mult <- function(mvpp_out, obs){
#   esout <- es_wrapper(mvpp_out, obs)
#   vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
#   vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
#   vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
#   vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
#   return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
# }


run_processing <- function(data, timeWindow, progress_ind = FALSE, compute_crps, ...){
  
  # Stations and days from the data
  stations <- unique(data$stat)
  days <- sort(unique(data$td))
  nout <- length(days)
  
  # 1 weather variable over 1 lead time with d stations
  d <- length(stations)
  
  # Ensemble members
  m <- sum(grepl("laef", names(data)))
  ensembleMembers <- sapply(1:m, FUN = function(x) paste0("laef", x))
  
  # generate objects to save scores to
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh", "gca","clayton","frank","gumbel")
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- param_list <- indep_list <- list()
  
  # Timing_list stores the time needed to do all relevant computations
  timing_list <- list()
  
  for(mm in 1:length(modelnames)){
    es_list[[mm]] <- vs1_list[[mm]] <- vs1w_list[[mm]] <- 
      vs0_list[[mm]] <- vs0w_list[[mm]] <- param_list[[mm]] <- indep_list[[mm]] <- matrix(NA, nrow = 1, ncol = nout) 
    crps_list[[mm]] <- array(NA, dim = c(nout, d))
    timing_list[[mm]] <- array(NA, dim = 1) 
  }
  names(crps_list) <- names(es_list) <- names(vs1_list) <- names(vs1w_list) <- 
    names(vs0_list) <- names(vs0w_list) <- names(param_list) <- names(indep_list) <- names(timing_list) <- modelnames
  


  # set random seed
  set.seed(1)
  
  # Create ensfc data structure
  ensfc <- array(NA, dim = c(nout, m, d))
  
  for (station in stations) {
    index <- match(station,stations)
    
    for (day in 1:length(days)) {
      # Extract forecast for day and dim = index
      try(
        ensfc[day,,index] <- unlist(subset(data, td == days[day] & stat == station)[ensembleMembers], use.names = FALSE)
      )
    }
  }
  
  # .GlobalEnv$test <- ensfc
  # postprocess ensemble forecasts
  start_time <- Sys.time()
  # emos_param[nn,dd,] <- c(loc, sc)
  pp_out <- array(NA, dim = c(length(days), d, 2))
  for (station in stations) {
    
    # Uses 50 day time window
    emos.result <- emos_T2M_mean_singleForecast(subset(data, stat == station))
    
    for (row in 1:length(emos.result$obs)) {
      values <- unlist(unname(emos.result[row,][c("ens_mu", "ens_sd")]))
      day <- match(emos.result[row,]$td, days)
      pp_out[day,match(station,stations),] <- values
    }
 

  }
  .GlobalEnv$test <- pp_out
  end_time <- Sys.time()
  
  return(pp_out)
  

}
k <- run_processing(data1, 3, TRUE, TRUE)





  # # if there are NaN's in pp output, generate new sets of obs and fc, and re-run pp code
  # # only happens very rarely for extreme parameter combinations
  # while(anyNA(pp_out)){
  #   set.seed(sample(1:100,1))
  #   obs <- generate_obs(model = obsmodel, nout = nout, ninit = ninit, ...)
  #   fc <- generate_ensfc(model = fcmodel, nout = nout, ninit = ninit, nmembers = nmembers, ...)
  #   pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init, 
  #                      obs = obs$obs, obs_init = obs$obs_init, 
  #                      train = "init", trainlength = NULL, emos_plus = TRUE)
  #   print("Repeating observations")
  # }
  # 
  # ## iterate over models and compute scores
  # 
  # # ensemble forecasts
  # start_time <- Sys.time()
  # if(compute_crps){
  #   crps_list$ens[rr, , ] <- crps_wrapper(fc$ensfc, obs$obs)
  # }
  # tmp <- eval_all_mult(mvpp_out = fc$ensfc, obs = obs$obs)
  # es_list$ens[rr, ] <- tmp$es 
  # vs1_list$ens[rr, ] <- tmp$vs1
  # vs1w_list$ens[rr, ] <- tmp$vs1w
  # vs0_list$ens[rr, ] <- tmp$vs0
  # vs0w_list$ens[rr, ] <- tmp$vs0w
  # 
  # end_time <- Sys.time()
  # 
  # timing_list$ens[rr] <- end_time - start_time
  # 
  # 
  # # EMOS.Q
  # start_time <- Sys.time()
  # emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #                obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow)
  # 
  # 
  # 
  # if(compute_crps){
  #   crps_list$emos.q[rr, , ] <- crps_wrapper(emos.q$mvppout, obs$obs)
  # }
  # 
  # tmp <- eval_all_mult(mvpp_out = emos.q$mvppout, obs = obs$obs)
  # es_list$emos.q[rr, ] <- tmp$es 
  # vs1_list$emos.q[rr, ] <- tmp$vs1
  # vs1w_list$emos.q[rr, ] <- tmp$vs1w
  # vs0_list$emos.q[rr, ] <- tmp$vs0
  # vs0w_list$emos.q[rr, ] <- tmp$vs0w
  # 
  # end_time <- Sys.time()
  # 
  # timing_list$emos.q[rr] <- end_time - start_time
  # 
  # # ECC.Q
  # start_time <- Sys.time()
  # ecc.q <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #               obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, 
  #               EMOS_sample = emos.q$mvppout, timeWindow = timeWindow)
  # 
  # if(compute_crps){
  #   crps_list$ecc.q[rr, , ] <- crps_wrapper(ecc.q$mvppout, obs$obs)
  # }
  # tmp <- eval_all_mult(mvpp_out = ecc.q$mvppout, obs = obs$obs)
  # es_list$ecc.q[rr, ] <- tmp$es 
  # vs1_list$ecc.q[rr, ] <- tmp$vs1
  # vs1w_list$ecc.q[rr, ] <- tmp$vs1w
  # vs0_list$ecc.q[rr, ] <- tmp$vs0
  # vs0w_list$ecc.q[rr, ] <- tmp$vs0w
  # 
  # end_time <- Sys.time()
  # 
  # timing_list$ecc.q[rr] <- end_time - start_time
  # 
  # # ECC.S -> involves randomness -> repeat rand_rep times
  # es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
  #   vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  # crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  # timing_list_tmp <- array(NA, dim = rand_rep)
  # for(RR in 1:rand_rep){
  #   start_time <- Sys.time()
  #   emos.s <- mvpp(method = "EMOS", variant = "S", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow)
  #   ecc.s <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #                 obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, 
  #                 EMOS_sample = emos.s$mvppout, timeWindow = timeWindow)
  #   
  #   if(compute_crps){
  #     crps_list_tmp[,,RR] <- crps_wrapper(ecc.s$mvppout, obs$obs)
  #   }
  #   tmp <- eval_all_mult(mvpp_out = ecc.s$mvppout, obs = obs$obs)
  #   es_list_tmp[,RR] <- tmp$es
  #   vs1_list_tmp[,RR] <- tmp$vs1
  #   vs1w_list_tmp[,RR] <- tmp$vs1w
  #   vs0_list_tmp[,RR] <- tmp$vs0
  #   vs0w_list_tmp[,RR] <- tmp$vs0w
  #   
  #   end_time <- Sys.time()
  #   
  #   timing_list_tmp[RR] <- end_time - start_time
  # }
  # crps_list$ecc.s[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
  # es_list$ecc.s[rr, ] <- apply(es_list_tmp, 1, mean)
  # vs1_list$ecc.s[rr, ] <- apply(vs1_list_tmp, 1, mean)
  # vs1w_list$ecc.s[rr, ] <- apply(vs1w_list_tmp, 1, mean)
  # vs0_list$ecc.s[rr, ] <- apply(vs0_list_tmp, 1, mean)
  # vs0w_list$ecc.s[rr, ] <- apply(vs0w_list_tmp, 1, mean)
  # timing_list$ecc.s[rr] <- apply(timing_list_tmp, 1, mean)
  # 
  # # dECC.Q
  # start_time <- Sys.time()
  # decc.q <- mvpp(method = "dECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #                obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, 
  #                EMOS_sample = emos.q$mvppout, ECC_out = ecc.q$mvppout, timeWindow = timeWindow) 
  # if(compute_crps){
  #   crps_list$decc.q[rr, , ] <- crps_wrapper(decc.q$mvppout, obs$obs)
  # }
  # tmp <- eval_all_mult(mvpp_out = decc.q$mvppout, obs = obs$obs)
  # es_list$decc.q[rr, ] <- tmp$es 
  # vs1_list$decc.q[rr, ] <- tmp$vs1
  # vs1w_list$decc.q[rr, ] <- tmp$vs1w
  # vs0_list$decc.q[rr, ] <- tmp$vs0
  # vs0w_list$decc.q[rr, ] <- tmp$vs0w
  # 
  # end_time <- Sys.time()
  # 
  # timing_list$decc.q[rr] <- end_time - start_time
  # 
  # # SSh -> involves randomness -> repeat rand_rep times
  # es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
  #   vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  # crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  # timing_list_tmp <- array(NA, dim = rand_rep)
  # for(RR in 1:rand_rep){
  #   start_time <- Sys.time()
  #   ssh <- mvpp(method = "SSh", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #               obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
  #               EMOS_sample = emos.q$mvppout, timeWindow = timeWindow) 
  #   
  #   if(compute_crps){
  #     crps_list_tmp[,,RR] <- crps_wrapper(ssh$mvppout, obs$obs)
  #   }
  #   tmp <- eval_all_mult(mvpp_out = ssh$mvppout, obs = obs$obs)
  #   es_list_tmp[,RR] <- tmp$es
  #   vs1_list_tmp[,RR] <- tmp$vs1
  #   vs1w_list_tmp[,RR] <- tmp$vs1w
  #   vs0_list_tmp[,RR] <- tmp$vs0
  #   vs0w_list_tmp[,RR] <- tmp$vs0w
  #   
  #   end_time <- Sys.time()
  #   
  #   timing_list_tmp[RR] <- end_time - start_time
  # }
  # crps_list$ssh[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
  # es_list$ssh[rr, ] <- apply(es_list_tmp, 1, mean)
  # vs1_list$ssh[rr, ] <- apply(vs1_list_tmp, 1, mean)
  # vs1w_list$ssh[rr, ] <- apply(vs1w_list_tmp, 1, mean)
  # vs0_list$ssh[rr, ] <- apply(vs0_list_tmp, 1, mean)
  # vs0w_list$ssh[rr, ] <- apply(vs0w_list_tmp, 1, mean)
  # timing_list$ssh[rr] <- apply(timing_list_tmp, 1, mean)
  # 
  # # GCA -> involves randomness -> repeat rand_rep times
  # es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
  #   vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  # crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  # timing_list_tmp <- array(NA, dim = rand_rep)
  # for(RR in 1:rand_rep){
  #   start_time <- Sys.time()
  #   gca <- mvpp(method = "GCA", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #               obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow) 
  #   if(compute_crps){
  #     crps_list_tmp[,,RR] <- crps_wrapper(gca$mvppout, obs$obs)
  #   }
  #   tmp <- eval_all_mult(mvpp_out = gca$mvppout, obs = obs$obs)
  #   es_list_tmp[,RR] <- tmp$es
  #   vs1_list_tmp[,RR] <- tmp$vs1
  #   vs1w_list_tmp[,RR] <- tmp$vs1w
  #   vs0_list_tmp[,RR] <- tmp$vs0
  #   vs0w_list_tmp[,RR] <- tmp$vs0w
  #   
  #   end_time <- Sys.time()
  #   
  #   timing_list_tmp[RR] <- end_time - start_time
  # }
  # crps_list$gca[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
  # es_list$gca[rr, ] <- apply(es_list_tmp, 1, mean)
  # vs1_list$gca[rr, ] <- apply(vs1_list_tmp, 1, mean)
  # vs1w_list$gca[rr, ] <- apply(vs1w_list_tmp, 1, mean)
  # vs0_list$gca[rr, ] <- apply(vs0_list_tmp, 1, mean)
  # vs0w_list$gca[rr, ] <- apply(vs0w_list_tmp, 1, mean)
  # timing_list$gca[rr] <- apply(timing_list_tmp, 1, mean)
  # 
  # # Archimedean copulas -> involves randomness -> repeat rand_rep times
  # for (method in c("Clayton","Frank", "Gumbel")) {
  #   es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
  #     vs0_list_tmp <- vs0w_list_tmp <- param_list_temp <- indep_list_temp <- matrix(NA, nrow = nout, ncol = rand_rep)
  #   crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  #   timing_list_tmp <- array(NA, dim = rand_rep)
  #   for(RR in 1:rand_rep){
  #     start_time <- Sys.time()
  #     mvd <- mvpp(method = method, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
  #                 obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, timeWindow = timeWindow) 
  #     
  #     if(compute_crps){
  #       crps_list_tmp[,,RR] <- crps_wrapper(mvd$mvppout, obs$obs)
  #     }
  #     tmp <- eval_all_mult(mvpp_out = mvd$mvppout, obs = obs$obs)
  #     es_list_tmp[,RR] <- tmp$es
  #     vs1_list_tmp[,RR] <- tmp$vs1
  #     vs1w_list_tmp[,RR] <- tmp$vs1w
  #     vs0_list_tmp[,RR] <- tmp$vs0
  #     vs0w_list_tmp[,RR] <- tmp$vs0w
  #     param_list_temp[,RR] <- mvd$params
  #     indep_list_temp[,RR] <- mvd$indep
  #     
  #     
  #     end_time <- Sys.time()
  #     
  #     timing_list_tmp[RR] <- end_time - start_time
  #   }
  #   
  #   if (method == "Clayton") {
  #     crps_list$clayton[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
  #     es_list$clayton[rr, ] <- apply(es_list_tmp, 1, mean)
  #     vs1_list$clayton[rr, ] <- apply(vs1_list_tmp, 1, mean)
  #     vs1w_list$clayton[rr, ] <- apply(vs1w_list_tmp, 1, mean)
  #     vs0_list$clayton[rr, ] <- apply(vs0_list_tmp, 1, mean)
  #     vs0w_list$clayton[rr, ] <- apply(vs0w_list_tmp, 1, mean)
  #     param_list$clayton[rr, ] <- apply(param_list_temp, 1, mean)
  #     timing_list$clayton[rr] <- apply(timing_list_tmp, 1, mean)
  #     indep_list$clayton[rr,] <- apply(indep_list_temp, 1, mean)
  #   } else if (method == "Frank") {
  #     crps_list$frank[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
  #     es_list$frank[rr, ] <- apply(es_list_tmp, 1, mean)
  #     vs1_list$frank[rr, ] <- apply(vs1_list_tmp, 1, mean)
  #     vs1w_list$frank[rr, ] <- apply(vs1w_list_tmp, 1, mean)
  #     vs0_list$frank[rr, ] <- apply(vs0_list_tmp, 1, mean)
  #     vs0w_list$frank[rr, ] <- apply(vs0w_list_tmp, 1, mean)
  #     param_list$frank[rr, ] <- apply(param_list_temp, 1, mean)
  #     timing_list$frank[rr] <- apply(timing_list_tmp, 1, mean)
  #     indep_list$frank[rr,] <- apply(indep_list_temp, 1, mean)
  #   } else if (method == "Gumbel") {
  #     crps_list$gumbel[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
  #     es_list$gumbel[rr, ] <- apply(es_list_tmp, 1, mean)
  #     vs1_list$gumbel[rr, ] <- apply(vs1_list_tmp, 1, mean)
  #     vs1w_list$gumbel[rr, ] <- apply(vs1w_list_tmp, 1, mean)
  #     vs0_list$gumbel[rr, ] <- apply(vs0_list_tmp, 1, mean)
  #     vs0w_list$gumbel[rr, ] <- apply(vs0w_list_tmp, 1, mean)
  #     param_list$gumbel[rr, ] <- apply(param_list_temp, 1, mean)
  #     timing_list$gumbel[rr] <- apply(timing_list_tmp, 1, mean)
  #     indep_list$gumbel[rr,] <- apply(indep_list_temp, 1, mean)
  #   } 
  # }
  #   
  # 
  # 
  # # return results, as a huge list
  # out <- list("crps_list" = crps_list, "es_list" = es_list, "vs1_list" = vs1_list,
  #             "vs1w_list" = vs1w_list, "vs0_list" = vs0_list, "vs0w_list" = vs0w_list, "param_list" = param_list, "indep_list" = indep_list, "timing_list" = timing_list)
  # return(out)



# run_wrapper <- function(runID){
#   # sink(file = paste0(Rout_dir, "_model_",modelSetting,"_ID_", runID, ".Rout"))
#   # Frank copula can only deal with positive tau values
#   tau <- runif(1, 0, 1)
#   par_values <- as.numeric(input_par[ID, ])
#   
#   res <- do.call("run_setting1",c(obsmodel = observationsModel, fcmodel = forecastModel, nout = evalDays, ninit = trainingDays, 
#                                   nmembers = ensembleMembers, timeWindow = timeWindow,
#                                   MCrep = MC_reps, rand_rep = randomRepetitions, 
#                                   progress_ind = TRUE, compute_crps = TRUE,
#                                   input_par[rownames(input_par) == runID,],
#                                   tau = tau))
#   if (modelSetting == 2) {
#     res$tau <- tau
#   }
#   
#   savename <- paste0(Rdata_dir,"_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,"_ID_", runID, ".Rdata")
#   save(res, input_par, file = savename)
#   # sink()
# }
# 
# 
# # Run the simulation for all parameter coefficients with a unique ID
# for (ID in 1:dim(input_par)[1]) {
#   closeAllConnections()
#   print(ID)
#   run_wrapper(runID = ID)
# }
# 


emos_T2M_mean_singleForecast(data1)


