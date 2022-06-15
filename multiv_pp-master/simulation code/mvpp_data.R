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
source("getData.R")
getData()






eval_all_mult <- function(mvpp_out, obs){
  esout <- es_wrapper(mvpp_out, obs)
  vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
  vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
  vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
  vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
  return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
}


run_processing <- function(data, trainingDays, progress_ind = FALSE, ...){
  
  # Use same time window as training days for UVPP
  timeWindow <- trainingDays
  
  # Stations and days from the data
  stations <- unique(data$stat)
  days <- sort(unique(data$td))
  nout <- length(days) - trainingDays
  
  # 1 weather variable over 1 lead time with d stations
  d <- length(stations)
  
  # Ensemble members
  m <- sum(grepl("laef", names(data)))
  ensembleMembers <- sapply(1:m, FUN = function(x) paste0("laef", x))
  
  # generate objects to save scores to
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh", "gca","clayton","frank","gumbel")
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- mvpp_list <- list()
  
  # Timing_list stores the time needed to do all relevant computations
  timing_list <- list()
  
  for(mm in 1:length(modelnames)){
    es_list[[mm]] <- vs1_list[[mm]] <- vs1w_list[[mm]] <- 
      vs0_list[[mm]] <- vs0w_list[[mm]] <- matrix(NA, nrow = 1, ncol = nout) 
    crps_list[[mm]] <- array(NA, dim = c(nout, d))
    timing_list[[mm]] <- array(NA, dim = 1) 
  }
  names(crps_list) <- names(es_list) <- names(vs1_list) <- names(vs1w_list) <- 
    names(vs0_list) <- names(vs0w_list) <- names(timing_list) <- modelnames
  


  # set random seed
  set.seed(1)
  
  # Create ensfc and obs data structure
  ensfc_init <- array(NA, dim = c(trainingDays, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  obs_init <- array(NA, dim = c(trainingDays, d))
  obs <- array(NA, dim = c(nout, d))
  
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
            ensfc_init[day,,index] <- unlist(dat[ensembleMembers], use.names = FALSE)
            obs_init[day, index] <- unlist(dat["obs"], use.names = FALSE)
          } else{
            ensfc[day - trainingDays,,index] <- unlist(dat[ensembleMembers], use.names = FALSE)
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
  
  
  
  print("Start UVPP")
  ## Start UVPP
  
  start_time <- Sys.time()
  uvpp <- list()
  
  for (station in stations) {
    
    statIndex <- match(station,stations)

    emos.result <- emos_T2M_mean_singleForecast(subset(data, stat == station), trainingDays)
    emos.result$stat <- statIndex

    uvpp <- rbind(uvpp, emos.result)
 

  }
  
  end_time <- Sys.time()
  
  timing_list$uvpp <- end_time - start_time
  
  print("Start MVPP")
  ## Start MVPP
  
  start_time <- Sys.time()
  
  dat <- uvpp[c("crps_emos", "stat")]
  
  # Create pp_out data structure and fill crps scores
  pp_out <- array(NA, dim = c(nout, d, 2))
  for (station in stations) {
    # Element of 1..d
    index <- match(station,stations)
    
    # Reformat days to be in 1..days
    for (day in (trainingDays+1):length(days)) {
      tryCatch(
        {
          temp <- subset(uvpp, td == days[day] & stat == index)
          crps_list$ens[day - trainingDays, index] <- unlist(unname(temp["crps_emos"]))
          pp_out[day - trainingDays, index, ] <- unlist(unname(temp[c("ens_mu", "ens_sd")]))
      }
      ,
      error = function(cond) {
        print(cond)
      })
    }
  }

  tmp <- eval_all_mult(mvpp_out = ensfc, obs = obs)
  es_list$ens <- tmp$es
  vs1_list$ens <- tmp$vs1
  vs1w_list$ens <- tmp$vs1w
  vs0_list$ens <- tmp$vs0
  vs0w_list$ens <- tmp$vs0w
  
  .GlobalEnv$uvpp <- uvpp
  .GlobalEnv$pp_out <- pp_out
  .GlobalEnv$obs_init <- obs_init
  .GlobalEnv$obs <- obs
  .GlobalEnv$ensfc_init <- ensfc_init
  .GlobalEnv$ensfc <- ensfc
  
  mvpp_list$ens <- ensfc
  
  print("EMOS.Q")
  # EMOS.Q
  start_time <- Sys.time()
  emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow)



  crps_list$emos.q <- crps_wrapper(emos.q$mvppout, obs)

  tmp <- eval_all_mult(mvpp_out = emos.q$mvppout, obs = obs)
  es_list$emos.q <- tmp$es
  vs1_list$emos.q <- tmp$vs1
  vs1w_list$emos.q <- tmp$vs1w
  vs0_list$emos.q <- tmp$vs0
  vs0w_list$emos.q <- tmp$vs0w

  end_time <- Sys.time()
  
  timing_list$emos.q <- end_time - start_time
  
  mvpp_list$emos.q <- emos.q$mvppout
  
  print("ECC.Q")
  # ECC.Q
  start_time <- Sys.time()
  ecc.q <- mvpp(method = "ECC", ensfc = ensfc, ensfc_init = ensfc_init,
                obs = obs, obs_init = obs_init, postproc_out = pp_out,
                EMOS_sample = emos.q$mvppout, timeWindow = timeWindow)

  crps_list$ecc.q <- crps_wrapper(ecc.q$mvppout, obs)
  tmp <- eval_all_mult(mvpp_out = ecc.q$mvppout, obs = obs)
  es_list$ecc.q <- tmp$es
  vs1_list$ecc.q <- tmp$vs1
  vs1w_list$ecc.q <- tmp$vs1w
  vs0_list$ecc.q <- tmp$vs0
  vs0w_list$ecc.q <- tmp$vs0w
  
  end_time <- Sys.time()
  
  timing_list$ecc.q <- end_time - start_time
  
  mvpp_list$ecc.q <- ecc.q$mvppout
  
  rand_rep <- 1

  print("ECC.S")
  # ECC.S -> involves randomness
  es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
    vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  timing_list_tmp <- array(NA, dim = rand_rep)
  for(RR in 1:rand_rep){
    start_time <- Sys.time()
    emos.s <- mvpp(method = "EMOS", variant = "S", ensfc = ensfc, ensfc_init = ensfc_init,
                   obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow)
    ecc.s <- mvpp(method = "ECC", ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = obs, obs_init = obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.s$mvppout, timeWindow = timeWindow)

    crps_list_tmp[,,RR] <- crps_wrapper(ecc.s$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = ecc.s$mvppout, obs = obs)
    es_list_tmp[,RR] <- tmp$es
    vs1_list_tmp[,RR] <- tmp$vs1
    vs1w_list_tmp[,RR] <- tmp$vs1w
    vs0_list_tmp[,RR] <- tmp$vs0
    vs0w_list_tmp[,RR] <- tmp$vs0w

    end_time <- Sys.time()

    timing_list_tmp[RR] <- end_time - start_time
  }
  crps_list$ecc.s <- apply(crps_list_tmp, c(1,2), mean)
  es_list$ecc.s <- apply(es_list_tmp, 1, mean)
  vs1_list$ecc.s <- apply(vs1_list_tmp, 1, mean)
  vs1w_list$ecc.s <- apply(vs1w_list_tmp, 1, mean)
  vs0_list$ecc.s <- apply(vs0_list_tmp, 1, mean)
  vs0w_list$ecc.s <- apply(vs0w_list_tmp, 1, mean)
  timing_list$ecc.s <- apply(timing_list_tmp, 1, mean)
  
  mvpp_list$ecc.s <- ecc.s$mvppout
  
  print("dECC.Q")
  # dECC.Q
  start_time <- Sys.time()
  decc.q <- mvpp(method = "dECC", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = obs, obs_init = obs_init, postproc_out = pp_out,
                 EMOS_sample = emos.q$mvppout, ECC_out = ecc.q$mvppout, timeWindow = timeWindow)

  crps_list$decc.q <- crps_wrapper(decc.q$mvppout, obs)

  tmp <- eval_all_mult(mvpp_out = decc.q$mvppout, obs = obs)
  es_list$decc.q <- tmp$es
  vs1_list$decc.q <- tmp$vs1
  vs1w_list$decc.q <- tmp$vs1w
  vs0_list$decc.q <- tmp$vs0
  vs0w_list$decc.q <- tmp$vs0w

  end_time <- Sys.time()

  timing_list$decc.q <- end_time - start_time
  
  mvpp_list$decc.q <- decc.q$mvppout
  
  print("SSh")
  # SSh -> involves randomness -> repeat rand_rep times
  es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
    vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  timing_list_tmp <- array(NA, dim = rand_rep)
  for(RR in 1:rand_rep){
    start_time <- Sys.time()
    ssh <- mvpp(method = "SSh", ensfc = ensfc, ensfc_init = ensfc_init,
                obs = obs, obs_init = obs_init, postproc_out = pp_out,
                EMOS_sample = emos.q$mvppout, timeWindow = timeWindow)

    crps_list_tmp[,,RR] <- crps_wrapper(ssh$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = ssh$mvppout, obs = obs)
    es_list_tmp[,RR] <- tmp$es
    vs1_list_tmp[,RR] <- tmp$vs1
    vs1w_list_tmp[,RR] <- tmp$vs1w
    vs0_list_tmp[,RR] <- tmp$vs0
    vs0w_list_tmp[,RR] <- tmp$vs0w

    end_time <- Sys.time()

    timing_list_tmp[RR] <- end_time - start_time
  }
  crps_list$ssh <- apply(crps_list_tmp, c(1,2), mean)
  es_list$ssh <- apply(es_list_tmp, 1, mean)
  vs1_list$ssh <- apply(vs1_list_tmp, 1, mean)
  vs1w_list$ssh <- apply(vs1w_list_tmp, 1, mean)
  vs0_list$ssh <- apply(vs0_list_tmp, 1, mean)
  vs0w_list$ssh <- apply(vs0w_list_tmp, 1, mean)
  timing_list$ssh <- apply(timing_list_tmp, 1, mean)
  
  mvpp_list$ssh <- ssh$mvppout

  print("GCA")
  # GCA -> involves randomness -> repeat rand_rep times
  es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
    vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  timing_list_tmp <- array(NA, dim = rand_rep)
  for(RR in 1:rand_rep){
    start_time <- Sys.time()
    gca <- mvpp(method = "GCA", ensfc = ensfc, ensfc_init = ensfc_init,
                obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, uvpp = uvpp)

    crps_list_tmp[,,RR] <- crps_wrapper(gca$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = gca$mvppout, obs = obs)
    es_list_tmp[,RR] <- tmp$es
    vs1_list_tmp[,RR] <- tmp$vs1
    vs1w_list_tmp[,RR] <- tmp$vs1w
    vs0_list_tmp[,RR] <- tmp$vs0
    vs0w_list_tmp[,RR] <- tmp$vs0w

    end_time <- Sys.time()

    timing_list_tmp[RR] <- end_time - start_time
  }
  crps_list$gca <- apply(crps_list_tmp, c(1,2), mean)
  es_list$gca <- apply(es_list_tmp, 1, mean)
  vs1_list$gca <- apply(vs1_list_tmp, 1, mean)
  vs1w_list$gca <- apply(vs1w_list_tmp, 1, mean)
  vs0_list$gca <- apply(vs0_list_tmp, 1, mean)
  vs0w_list$gca <- apply(vs0w_list_tmp, 1, mean)
  timing_list$gca <- apply(timing_list_tmp, 1, mean)
  
  mvpp_list$gca <- gca$mvppout

  # Archimedean copulas -> involves randomness -> repeat rand_rep times
  for (method in c("Clayton","Frank", "Gumbel")) {
    print(method)
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
      vs0_list_tmp <- vs0w_list_tmp  <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    timing_list_tmp <- array(NA, dim = rand_rep)
    for(RR in 1:rand_rep){
      start_time <- Sys.time()
      mvd <- mvpp(method = method, ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, uvpp = uvpp)
      
      .GlobalEnv$mvd <- mvd
      
      good_rows <- !unname(apply(mvd$mvppout, MARGIN = 1, FUN = function(x) any(is.na(x))))
      score_mvppout <- mvd$mvppout[good_rows,,]
      score_obs <- obs[good_rows,]

      crps_list_tmp[good_rows,,RR] <- crps_wrapper(score_mvppout, score_obs)
      tmp <- eval_all_mult(mvpp_out = score_mvppout, obs = score_obs)
      es_list_tmp[good_rows,RR] <- tmp$es
      vs1_list_tmp[good_rows,RR] <- tmp$vs1
      vs1w_list_tmp[good_rows,RR] <- tmp$vs1w
      vs0_list_tmp[good_rows,RR] <- tmp$vs0
      vs0w_list_tmp[good_rows,RR] <- tmp$vs0w


      end_time <- Sys.time()

      timing_list_tmp[RR] <- end_time - start_time
    }

    if (method == "Clayton") {
      crps_list$clayton <- apply(crps_list_tmp, c(1,2), mean)
      es_list$clayton <- apply(es_list_tmp, 1, mean)
      vs1_list$clayton <- apply(vs1_list_tmp, 1, mean)
      vs1w_list$clayton <- apply(vs1w_list_tmp, 1, mean)
      vs0_list$clayton <- apply(vs0_list_tmp, 1, mean)
      vs0w_list$clayton <- apply(vs0w_list_tmp, 1, mean)
      timing_list$clayton <- apply(timing_list_tmp, 1, mean)
      mvpp_list$clayton <- mvd
    } else if (method == "Frank") {
      crps_list$frank <- apply(crps_list_tmp, c(1,2), mean)
      es_list$frank <- apply(es_list_tmp, 1, mean)
      vs1_list$frank <- apply(vs1_list_tmp, 1, mean)
      vs1w_list$frank <- apply(vs1w_list_tmp, 1, mean)
      vs0_list$frank <- apply(vs0_list_tmp, 1, mean)
      vs0w_list$frank <- apply(vs0w_list_tmp, 1, mean)
      timing_list$frank <- apply(timing_list_tmp, 1, mean)
      mvpp_list$frank <- mvd
    } else if (method == "Gumbel") {
      crps_list$gumbel <- apply(crps_list_tmp, c(1,2), mean)
      es_list$gumbel <- apply(es_list_tmp, 1, mean)
      vs1_list$gumbel <- apply(vs1_list_tmp, 1, mean)
      vs1w_list$gumbel <- apply(vs1w_list_tmp, 1, mean)
      vs0_list$gumbel <- apply(vs0_list_tmp, 1, mean)
      vs0w_list$gumbel <- apply(vs0w_list_tmp, 1, mean)
      timing_list$gumbel <- apply(timing_list_tmp, 1, mean)
      mvpp_list$gumbel <- mvd
    }
  }
  
  
  print("Returning Output")
  # return results, as a huge list
  out <- list("crps_list" = crps_list, "es_list" = es_list, "vs1_list" = vs1_list,
              "vs1w_list" = vs1w_list, "vs0_list" = vs0_list, "vs0w_list" = vs0w_list, 
              "timing_list" = timing_list, "mvpp_list" = mvpp_list, "obs" = obs)
  return(out)
}
traininDays <- 30
saveDir <- "./../Data/Rdata_LAEF/"

res <- run_processing(data1, traininDays, progress_ind = TRUE)
savename <- paste0(saveDir, "Res_group_1", ".Rdata")
save(res, file = savename)

res <- run_processing(data2, traininDays, progress_ind = TRUE)
savename <- paste0(saveDir, "Res_group_2", ".Rdata")
save(res, file = savename)

res <- run_processing(data3, traininDays, progress_ind = TRUE)
savename <- paste0(saveDir, "Res_group_3", ".Rdata")
save(res, file = savename)



