# rm(list=ls())

library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)
library("here")

here2 <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  if ("RStudio" %in% args) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    file_arg <- "--file="
    filepath <- sub(file_arg, "", grep(file_arg, args, value = TRUE))
    dirname(filepath)
  }
}


setwd(here2())

# "source" files for functions
dir <- "./sourceArchimedean/"
source(paste0(dir, "postprocess_ensfc_arch.R"))
source(paste0(dir, "mvpp_arch.R"))
source(paste0(dir, "evaluation_functions_arch.R"))
source(paste0(dir, "CopulaParameter.R"))
source(paste0(dir, "mvpp_data_copula_study.R"))
source(paste0(dir, "generate_trainingdays.R"))
source(paste0(dir, "mvpp_arch_fixed_training.R"))
source("ECC_T2M_Emos_subfunctions.R")


# Load data
source("getData.R")
if (!("data1" %in% names(.GlobalEnv))) {
  getData()
}

eval_all_mult <- function(mvpp_out, obs){
  esout <- es_wrapper(mvpp_out, obs)
  vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
  vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
  vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
  vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
  return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
}


run_processing <- function(data, trainingDays, progress_ind = FALSE, timeWindow, fix_training_days, training_days_method){
  
  # Stations and days from the data
  stations <- unique(data$stat)
  days <- sort(unique(data$td))
  nout <- length(days) - trainingDays
  
  # 1 weather variable over 1 lead time with d stations
  d <- length(stations)
  
  # Ensemble members
  m <- sum(grepl("laef", names(data)))
  ecc_m <- timeWindow
  ensembleMembers <- sapply(1:m, FUN = function(x) paste0("laef", x))
  
  if (fix_training_days) {
    training_df <- getTrainingDays(trainingDays, nout, ecc_m, training_days_method)
  }
  
  
  # generate objects to save scores to
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh.h", "ssh.i", "gca", "gca.cop", "Clayton","Frank","Gumbel", "Surv_Gumbel")
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- mvpp_list <- chosenCopula_list <- list()
  
  # Timing_list stores the time needed to do all relevant computations
  timing_list <- list()
  
  for(mm in 1:length(modelnames)){
    es_list[[mm]] <- vs1_list[[mm]] <- vs1w_list[[mm]] <- 
      vs0_list[[mm]] <- vs0w_list[[mm]] <- matrix(NA, nrow = 1, ncol = nout) 
    crps_list[[mm]] <- array(NA, dim = c(nout, d))
    timing_list[[mm]] <- chosenCopula_list[[mm]] <- array(NA, dim = 1) 
  }
  names(crps_list) <- names(es_list) <- names(vs1_list) <- names(vs1w_list) <- 
    names(vs0_list) <- names(vs0w_list) <- names(timing_list) <- names(chosenCopula_list) <- modelnames
  


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
  # Start UVPP

  start_time <- Sys.time()
  uvpp <- list()

  for (station in stations) {

    statIndex <- match(station,stations)

    emos.result <- emos_T2M_mean_singleForecast(subset(data, stat == station), trainingDays)
    emos.result$stat <- statIndex

    uvpp <- rbind(uvpp, emos.result)


  }

  .GlobalEnv$uvpp <- uvpp

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
  
  .GlobalEnv$uvpp <- uvpp
  .GlobalEnv$pp_out <- pp_out
  .GlobalEnv$obs_init <- obs_init
  .GlobalEnv$obs <- obs
  .GlobalEnv$ensfc_init <- ensfc_init
  .GlobalEnv$ensfc <- ensfc

  tmp <- eval_all_mult(mvpp_out = ensfc, obs = obs)
  es_list$ens <- tmp$es
  vs1_list$ens <- tmp$vs1
  vs1w_list$ens <- tmp$vs1w
  vs0_list$ens <- tmp$vs0
  vs0w_list$ens <- tmp$vs0w
  
  
  
  mvpp_list$ens <- ensfc
  
  print("EMOS.Q")
  # EMOS.Q
  start_time <- Sys.time()
  emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = m)



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
                   obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = m)
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
  
  
  
  print("SSh-H")
  # SSh -> involves randomness -> repeat rand_rep times
  es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
    vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  timing_list_tmp <- array(NA, dim = rand_rep)
  
  print("EMOS.Q with different m")
  # EMOS.Q with different m
  emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m)
  
  for(RR in 1:rand_rep){
    print(RR)
    start_time <- Sys.time()
    
    if (fix_training_days) {
      ssh.h <- mvpp_fixed(method = "SSh-H", ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = obs, obs_init = obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m, training_df = training_df)
    } else {
      ssh.h <- mvpp(method = "SSh-H", ensfc = ensfc, ensfc_init = ensfc_init,
                    obs = obs, obs_init = obs_init, postproc_out = pp_out,
                    EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m)
    }

    crps_list_tmp[,,RR] <- crps_wrapper(ssh.h$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = ssh.h$mvppout, obs = obs)
    es_list_tmp[,RR] <- tmp$es
    vs1_list_tmp[,RR] <- tmp$vs1
    vs1w_list_tmp[,RR] <- tmp$vs1w
    vs0_list_tmp[,RR] <- tmp$vs0
    vs0w_list_tmp[,RR] <- tmp$vs0w

    end_time <- Sys.time()

    timing_list_tmp[RR] <- end_time - start_time
  }
  crps_list$ssh.h <- apply(crps_list_tmp, c(1,2), mean)
  es_list$ssh.h <- apply(es_list_tmp, 1, mean)
  vs1_list$ssh.h <- apply(vs1_list_tmp, 1, mean)
  vs1w_list$ssh.h <- apply(vs1w_list_tmp, 1, mean)
  vs0_list$ssh.h <- apply(vs0_list_tmp, 1, mean)
  vs0w_list$ssh.h <- apply(vs0w_list_tmp, 1, mean)
  timing_list$ssh.h <- apply(timing_list_tmp, 1, mean)
  
  mvpp_list$ssh.h <- ssh.h$mvppout
  
  print("SSH-I14")
  
  for(RR in 1:rand_rep){
    start_time <- Sys.time()
    
    if (fix_training_days) {
      ssh.i <- mvpp_fixed(method = "SSh-I14", ensfc = ensfc, ensfc_init = ensfc_init,
                    obs = obs, obs_init = obs_init, postproc_out = pp_out,
                    EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m, training_df = training_df)
    } else {
      ssh.i <- mvpp(method = "SSh-I14", ensfc = ensfc, ensfc_init = ensfc_init,
                    obs = obs, obs_init = obs_init, postproc_out = pp_out,
                    EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m)
    }
    
    crps_list_tmp[,,RR] <- crps_wrapper(ssh.i$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = ssh.i$mvppout, obs = obs)
    es_list_tmp[,RR] <- tmp$es
    vs1_list_tmp[,RR] <- tmp$vs1
    vs1w_list_tmp[,RR] <- tmp$vs1w
    vs0_list_tmp[,RR] <- tmp$vs0
    vs0w_list_tmp[,RR] <- tmp$vs0w
    
    end_time <- Sys.time()
    
    timing_list_tmp[RR] <- end_time - start_time
  }
  crps_list$ssh.i <- apply(crps_list_tmp, c(1,2), mean)
  es_list$ssh.i <- apply(es_list_tmp, 1, mean)
  vs1_list$ssh.i <- apply(vs1_list_tmp, 1, mean)
  vs1w_list$ssh.i <- apply(vs1w_list_tmp, 1, mean)
  vs0_list$ssh.i <- apply(vs0_list_tmp, 1, mean)
  vs0w_list$ssh.i <- apply(vs0w_list_tmp, 1, mean)
  timing_list$ssh.i <- apply(timing_list_tmp, 1, mean)
  
  mvpp_list$ssh.i <- ssh.i$mvppout

  print("GCA")
  # GCA -> involves randomness -> repeat rand_rep times
  es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
    vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  timing_list_tmp <- array(NA, dim = rand_rep)
  for(RR in 1:rand_rep){
    start_time <- Sys.time()
    
    if (fix_training_days) {
      gca <- mvpp_fixed(method = "GCA", ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp, training_df = training_df)
    } else {
      gca <- mvpp(method = "GCA", ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
    }

    .GlobalEnv$tt <- gca$mvppout
    
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
  
  print("CopGCA")
  # GCA -> involves randomness -> repeat rand_rep times
  es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <-
    vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
  crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
  timing_list_tmp <- array(NA, dim = rand_rep)
  for(RR in 1:rand_rep){
    start_time <- Sys.time()
    
    if (fix_training_days) {
      gca.cop <- mvpp_fixed(method = "CopGCA", ensfc = ensfc, ensfc_init = ensfc_init,
                      obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp, training_df = training_df)
    } else {
      gca.cop <- mvpp(method = "CopGCA", ensfc = ensfc, ensfc_init = ensfc_init,
                      obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
    }
    
    .GlobalEnv$tt <- gca.cop$mvppout
    
    crps_list_tmp[,,RR] <- crps_wrapper(gca.cop$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = gca.cop$mvppout, obs = obs)
    es_list_tmp[,RR] <- tmp$es
    vs1_list_tmp[,RR] <- tmp$vs1
    vs1w_list_tmp[,RR] <- tmp$vs1w
    vs0_list_tmp[,RR] <- tmp$vs0
    vs0w_list_tmp[,RR] <- tmp$vs0w
    
    end_time <- Sys.time()
    
    timing_list_tmp[RR] <- end_time - start_time
  }
  crps_list$gca.cop <- apply(crps_list_tmp, c(1,2), mean)
  es_list$gca.cop <- apply(es_list_tmp, 1, mean)
  vs1_list$gca.cop <- apply(vs1_list_tmp, 1, mean)
  vs1w_list$gca.cop <- apply(vs1w_list_tmp, 1, mean)
  vs0_list$gca.cop <- apply(vs0_list_tmp, 1, mean)
  vs0w_list$gca.cop <- apply(vs0w_list_tmp, 1, mean)
  timing_list$gca.cop <- apply(timing_list_tmp, 1, mean)
  
  mvpp_list$gca.cop <- gca.cop$mvppout

  # Archimedean copulas -> involves randomness -> repeat rand_rep times
  for (method in c("Surv_Gumbel", "Clayton","Frank", "Gumbel")) {
    print(method)


      
    start_time <- Sys.time()
    
    if (fix_training_days) {
      mvd <- mvpp_fixed(method = method, ensfc = ensfc, ensfc_init = ensfc_init,
                obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp, training_df = training_df)
    } else {
      mvd <- mvpp(method = method, ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
    }
    
    chosenCopula_list[[method]] <- mvd$chosenCopula
    
    # good_rows <- !unname(apply(mvd$mvppout, MARGIN = 1, FUN = function(x) any(is.na(x))))
    # score_mvppout <- mvd$mvppout[good_rows,,]
    # score_obs <- obs[good_rows,]

    crps_list[[method]] <- crps_wrapper(mvd$mvppout, obs)
    tmp <- eval_all_mult(mvpp_out = mvd$mvppout, obs = obs)
    es_list[[method]] <- tmp$es
    vs1_list[[method]] <- tmp$vs1
    vs1w_list[[method]] <- tmp$vs1w
    vs0_list[[method]] <- tmp$vs0
    vs0w_list[[method]] <- tmp$vs0w
    
    mvpp_list[[method]] <- mvd$mvppout

    end_time <- Sys.time()

    timing_list[[method]] <- end_time - start_time

  }
  
  # Seasonal copula MVPP
  # print("Seasonal")
  # start_time <- Sys.time()
  # 
  # mvd <- seasonal_copula(data, ensfc = ensfc, obs = obs, obs_init = obs_init, postproc_out = pp_out, ecc_m = ecc_m)
  # 
  # good_rows <- !unname(apply(mvd$mvppout, MARGIN = 1, FUN = function(x) any(is.na(x))))
  # score_mvppout <- mvd$mvppout[good_rows,,]
  # score_obs <- obs[good_rows,]
  # 
  # crps_list[["Seasonal"]][good_rows,] <- crps_wrapper(score_mvppout, score_obs)
  # tmp <- eval_all_mult(mvpp_out = score_mvppout, obs = score_obs)
  # es_list[["Seasonal"]][good_rows] <- tmp$es
  # vs1_list[["Seasonal"]][good_rows] <- tmp$vs1
  # vs1w_list[["Seasonal"]][good_rows] <- tmp$vs1w
  # vs0_list[["Seasonal"]][good_rows] <- tmp$vs0
  # vs0w_list[["Seasonal"]][good_rows] <- tmp$vs0w
  # 
  # mvpp_list[["Seasonal"]] <- mvd$mvppout
  # 
  # end_time <- Sys.time()
  # 
  # timing_list[["Seasonal"]] <- end_time - start_time
  # 
  # copula_df <- mvd$copula_df
  
  
  print("Returning Output")
  # return results, as a huge list
  out <- list("crps_list" = crps_list, "es_list" = es_list, "vs1_list" = vs1_list,
              "vs1w_list" = vs1w_list, "vs0_list" = vs0_list, "vs0w_list" = vs0w_list, 
              "timing_list" = timing_list, "mvpp_list" = mvpp_list, "chosenCopula_list" = chosenCopula_list,
              "obs" = obs)
  return(out)
}

# fix_training_days <- FALSE
# training_days_method <- "random_2w_interval"



compute_res <-function(trainingDays, timeWindow, fix_training_days, training_days_method) {
  saveDir <- "./../Data/Rdata_LAEF/"
  saveDir <- paste0(saveDir, "m_", timeWindow, "_")
  
  if (fix_training_days) {
    saveDir <- paste0(saveDir, training_days_method, "_")
  }
  
  print(paste0("Start processing for directory: ", saveDir))
  
  print("")
  print("Group 1")
  print("")
  res <- run_processing(data1, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_1", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 2")
  print("")
  res <- run_processing(data2, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_2", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 3")
  print("")
  res <- run_processing(data3, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_3", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 4")
  print("")
  res <- run_processing(data4, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_4", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 5")
  print("")
  res <- run_processing(data5, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_5", ".Rdata")
  save(res, file = savename)
}

trainingDays <- 365

for (fix_training_days in c(TRUE, FALSE)) {
  if (fix_training_days) {
    for (training_days_method in c("random_2w_interval")) {

      compute_res(trainingDays, 100, fix_training_days, training_days_method)

    }
  } else {
    compute_res(trainingDays, 100, fix_training_days, training_days_method)
  }
}


compute_res(trainingDays, 100, FALSE, training_days_method)
