# rm(list=ls())

library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)

# "source" files for functions
dir <- "simulation code/sourceArchimedean/"
source(paste0(dir, "postprocess_ensfc_arch.R"))
source(paste0(dir, "mvpp_arch.R"))
source(paste0(dir, "evaluation_functions_arch.R"))
source(paste0(dir, "CopulaParameter.R"))
source(paste0(dir, "mvpp_data_copula_study.R"))
source(paste0(dir, "generate_trainingdays.R"))
source(paste0(dir, "mvpp_arch_fixed_training.R"))
source("simulation code/ECC_T2M_Emos_subfunctions.R")


# Load data
for (groupNr in 1:5) {
  x <- load(paste0("Data/Input/data", groupNr, ".Rdata"))
  assign(paste0("data", groupNr), get(x))
  rm(list = x, envir = .GlobalEnv)
  
  x <- load(paste0("Data/UVPP/uvpp", groupNr, ".Rdata"))
  assign(paste0("uvpp", groupNr), get(x))
  rm(list = x, envir = .GlobalEnv)
  
  x <- load(paste0("Data/SimilarityMatrix/simMatrix", groupNr, ".Rdata"))
  assign(paste0("simMatrix", groupNr), get(x))
  rm(list = x, envir = .GlobalEnv)
}

eval_all_mult <- function(mvpp_out, obs) {
  esout <- es_wrapper(mvpp_out, obs)
  vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
  vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
  vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
  vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
  return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
}

run_processing <- function(data, uvpp, sim_matrix, trainingDays, progress_ind = FALSE, timeWindow, fix_training_days, training_days_method) {
  
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
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh.h", "ssh.i", "sim.ssh",
                  "gca", "gca.sh", "gca.cop", "gca.cop.sh", "Clayton", "Frank", "Gumbel", "Surv_Gumbel",
                  "Claytonsh", "Franksh", "Gumbelsh", "Surv_Gumbelsh")
  
  score.env <- new.env()
  score.env$crps_list <- score.env$es_list <- score.env$vs1_list <- score.env$vs1w_list <- score.env$vs0_list <- score.env$vs0w_list <- score.env$mvpp_list <- chosenCopula_list <- list()
  

  
  # Timing_list stores the time needed to do all relevant computations
  score.env$timing_list <- list()
  
  for (mm in 1:length(modelnames)){
    score.env$es_list[[mm]] <- score.env$vs1_list[[mm]] <- score.env$vs1w_list[[mm]] <- 
      score.env$vs0_list[[mm]] <- score.env$vs0w_list[[mm]] <- matrix(NA, nrow = 1, ncol = nout) 
    score.env$crps_list[[mm]] <- array(NA, dim = c(nout, d))
    score.env$timing_list[[mm]] <- chosenCopula_list[[mm]] <- array(NA, dim = 1) 
  }
  names(score.env$crps_list) <- names(score.env$es_list) <- names(score.env$vs1_list) <- names(score.env$vs1w_list) <- 
    names(score.env$vs0_list) <- names(score.env$vs0w_list) <- names(score.env$timing_list) <- names(chosenCopula_list) <- modelnames
  
  # set random seed
  set.seed(1)
  
  # Create ensfc and obs data structure
  ensfc_init <- array(NA, dim = c(trainingDays, m, d))
  ensfc <- array(NA, dim = c(nout, m, d))
  obs_init <- array(NA, dim = c(trainingDays, d))
  score.env$obs <- array(NA, dim = c(nout, d))
  
  # Dimension only consists of stations
  for (station in stations) {
    # Element of 1..d
    index <- match(station, stations)
    
    # Reformat days to be in 1..days
    for (day in 1:length(days)) {
      # Extract forecast for day and dim = index
      tryCatch(
        {
          dat <- subset(data, td == days[day] & stat == station)
          if (day <= trainingDays) {
            ensfc_init[day, , index] <- unlist(dat[ensembleMembers], use.names = FALSE)
            obs_init[day, index] <- unlist(dat["obs"], use.names = FALSE)
          } else {
            ensfc[day - trainingDays, , index] <- unlist(dat[ensembleMembers], use.names = FALSE)
            score.env$obs[day - trainingDays, index] <- unlist(dat["obs"], use.names = FALSE)
          }
          
        }
        ,
        error = function(cond) {
          print(cond)
        }
      )
    }
  }

  ##############
  ## Ensemble ##
  ##############

  print("Creating ensemble structure")
  
  start_time <- Sys.time()
  
  dat <- uvpp[c("crps_emos", "stat")]
  
  # Create pp_out data structure and fill crps scores
  pp_out <- array(NA, dim = c(nout, d, 2))
  for (station in stations) {
    # Element of 1..d
    index <- match(station, stations)
    
    # Reformat days to be in 1..days
    for (day in (trainingDays + 1):length(days)) {
      tryCatch(
               {
                 temp <- subset(uvpp, td == days[day] & stat == index)
                 score.env$crps_list$ens[day - trainingDays, index] <- unlist(unname(temp["crps_emos"]))
                 pp_out[day - trainingDays, index, ] <- unlist(unname(temp[c("ens_mu", "ens_sd")]))
               }
               ,
               error = function(cond) {
                 print(cond)
               })
    }
  }

  ########## 
  ## MVPP ##  
  ########## 

  print("Start MVPP")
  
  add_scores <- function(mvpp_out, method, start_time) {
    
    # CRPS
    score.env$crps_list[[method]] <- crps_wrapper(mvpp_out, score.env$obs)
    
    # ES and VS
    tmp <- eval_all_mult(mvpp_out = mvpp_out, obs = score.env$obs)
    score.env$es_list[[method]]   <- tmp$es
    score.env$vs1_list[[method]]  <- tmp$vs1
    score.env$vs1w_list[[method]] <- tmp$vs1w
    score.env$vs0_list[[method]]  <- tmp$vs0
    score.env$vs0w_list[[method]] <- tmp$vs0w
    score.env$mvpp_list[[method]] <- mvpp_out

    # Timing
    end_time <- Sys.time()
    score.env$timing_list[[method]] <- end_time - start_time
  }

  add_scores(ensfc, "ens", start_time)

  ############
  ## EMOS-Q ##
  ############
  
  print("EMOS.Q")

  start_time <- Sys.time()
  
  emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = m)

  add_scores(emos.q$mvppout, "emos.q", start_time)
  
  ###########
  ## ECC-Q ##
  ###########
  
  print("ECC.Q")

  start_time <- Sys.time()

  ecc.q <- mvpp(method = "ECC", ensfc = ensfc, ensfc_init = ensfc_init,
                obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                EMOS_sample = emos.q$mvppout, timeWindow = timeWindow)

  add_scores(ecc.q$mvppout, "ecc.q", start_time)
  
  ###########
  ## ECC-S ##
  ###########

  print("ECC.S")

  start_time <- Sys.time()

  emos.s <- mvpp(method = "EMOS", variant = "S", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = m)
  ecc.s <- mvpp(method = "ECC", ensfc = ensfc, ensfc_init = ensfc_init,
                obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                EMOS_sample = emos.s$mvppout, timeWindow = timeWindow)
  
  add_scores(ecc.s$mvppout, "ecc.s", start_time)

  ############
  ## dECC-Q ##
  ############
  
  print("dECC.Q")
  
  start_time <- Sys.time()

  decc.q <- mvpp(method = "dECC", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                 EMOS_sample = emos.q$mvppout, ECC_out = ecc.q$mvppout, timeWindow = timeWindow)

  add_scores(decc.q$mvppout, "decc.q", start_time)

  ###########
  ## SSH-H ##
  ###########
  
  print("SSh-H")
  
  print("EMOS.Q with different m")

  emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m)
  
    
  if (fix_training_days) {
    ssh.h <- mvpp_fixed(method = "SSh-H", ensfc = ensfc, ensfc_init = ensfc_init,
                        obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                        EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m, training_df = training_df)
  } else {
    ssh.h <- mvpp(method = "SSh-H", ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m)
  }

  add_scores(ssh.h$mvppout, "ssh.h", start_time)

  ###########
  ## SSH-I ##
  ###########

  print("SSH-I14")
  
  start_time <- Sys.time()
  
  if (fix_training_days) {
    ssh.i <- mvpp_fixed(method = "SSh-I14", ensfc = ensfc, ensfc_init = ensfc_init,
                        obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                        EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m, training_df = training_df)
  } else {
    ssh.i <- mvpp(method = "SSh-I14", ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m)
  }
    
  add_scores(ssh.i$mvppout, "ssh.i", start_time)
  
  ################
  ## SimSchaake ##
  ################
  
  print("SimSchaake")
  
  start_time <- Sys.time()
  
  sim.ssh <- mvpp(method = "SimSchaake", ensfc = ensfc, ensfc_init = ensfc_init, sim_matrix = sim_matrix,
                  obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, uvpp = uvpp,
                  EMOS_sample = emos.q$mvppout, timeWindow = timeWindow, ecc_m = ecc_m)
  
  .GlobalEnv$sim.ssh <- sim.ssh

  add_scores(sim.ssh$mvppout, "sim.ssh", start_time)

  #########
  ## GCA ##
  #########

  print("GCA")

  start_time <- Sys.time()
  
  if (fix_training_days) {
    gca <- mvpp_fixed(method = "GCA", ensfc = ensfc, ensfc_init = ensfc_init,
                      obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp, training_df = training_df)
  } else {
    gca <- mvpp(method = "GCA", ensfc = ensfc, ensfc_init = ensfc_init,
                obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
  }

  add_scores(gca$mvppout, "gca", start_time)

  ###########
  ## GCAsh ##
  ###########

  print("GCAsh")

  start_time <- Sys.time()

  gca.sh <- mvpp(method = "GCAsh", ensfc = ensfc, ensfc_init = ensfc_init,
                 obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, EMOS_sample = emos.q$mvppout, 
                 timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
    
  add_scores(gca.sh$mvppout, "gca.sh", start_time)

  ############
  ## CopGCA ##
  ############
  
  print("CopGCA")

  start_time <- Sys.time()
  
  if (fix_training_days) {
    gca.cop <- mvpp_fixed(method = "CopGCA", ensfc = ensfc, ensfc_init = ensfc_init,
                          obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp, training_df = training_df)
  } else {
    gca.cop <- mvpp(method = "CopGCA", ensfc = ensfc, ensfc_init = ensfc_init,
                    obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
  }
    
  add_scores(gca.cop$mvppout, "gca.cop", start_time)

  ##############
  ## CopGCAsh ##
  ##############
  
  print("CopGCAsh")

  start_time <- Sys.time()
  
  gca.cop.sh <- mvpp(method = "CopGCAsh", ensfc = ensfc, ensfc_init = ensfc_init,
                     obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, EMOS_sample = emos.q$mvppout, 
                     timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
    
  add_scores(gca.cop.sh$mvppout, "gca.cop.sh", start_time)

  #########################
  ## Archimedean Copulas ##
  #########################

  for (method in c("Surv_Gumbel", "Clayton", "Frank", "Gumbel", "Surv_Gumbelsh", "Claytonsh", "Franksh", "Gumbelsh")) {
    print(method)

    start_time <- Sys.time()
    
    if (fix_training_days) {
      mvd <- mvpp_fixed(method = method, ensfc = ensfc, ensfc_init = ensfc_init,
                        obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, timeWindow = timeWindow, 
                        ecc_m = ecc_m, uvpp = uvpp, training_df = training_df)
    } else {
      mvd <- mvpp(method = method, ensfc = ensfc, ensfc_init = ensfc_init,
                  obs = score.env$obs, obs_init = obs_init, postproc_out = pp_out, EMOS_sample = emos.q$mvppout, 
                  timeWindow = timeWindow, ecc_m = ecc_m, uvpp = uvpp)
    }
    
    chosenCopula_list[[method]] <- mvd$chosenCopula
    
    add_scores(mvd$mvppout, method, start_time)

  }

  print("Returning Output")
  # return results, as a huge list
  out <- list("crps_list" = score.env$crps_list, "es_list" = score.env$es_list, "vs1_list" = score.env$vs1_list,
              "vs1w_list" = score.env$vs1w_list, "vs0_list" = score.env$vs0_list, "vs0w_list" = score.env$vs0w_list, 
              "timing_list" = score.env$timing_list, "mvpp_list" = score.env$mvpp_list, "chosenCopula_list" = chosenCopula_list,
              "obs" = score.env$obs)
  return(out)
}

# Computes the results for all groups
compute_res <- function(trainingDays, timeWindow, fix_training_days, training_days_method) {
  saveDir <- "Data/Rdata_LAEF/"
  saveDir <- paste0(saveDir, "m_", timeWindow, "_")
  
  if (fix_training_days) {
    saveDir <- paste0(saveDir, training_days_method, "_")
  }
  
  print(paste0("Start processing for directory: ", saveDir))
  
  # print("")
  # print("Group 1")
  # print("")
  # res <- run_processing(data1, uvpp1, simMatrix1, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  # savename <- paste0(saveDir, "Res_group_1", ".Rdata")
  # save(res, file = savename)
  
  print("")
  print("Group 2")
  print("")
  res <- run_processing(data2, uvpp2, simMatrix2, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_2", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 3")
  print("")
  res <- run_processing(data3, uvpp3, simMatrix3, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_3", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 4")
  print("")
  res <- run_processing(data4, uvpp4, simMatrix4, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_4", ".Rdata")
  save(res, file = savename)
  
  print("")
  print("Group 5")
  print("")
  res <- run_processing(data5, uvpp5, simMatrix5, trainingDays, progress_ind = TRUE, timeWindow, fix_training_days, training_days_method)
  savename <- paste0(saveDir, "Res_group_5", ".Rdata")
  save(res, file = savename)
}

trainingDays <- 365

# for (fix_training_days in c(TRUE, FALSE)) {
#   if (fix_training_days) {
#     for (training_days_method in c("random_2w_interval")) {
# 
#       compute_res(trainingDays, 100, fix_training_days, training_days_method)
# 
#     }
#   } else {
#     compute_res(trainingDays, 100, fix_training_days, training_days_method)
#   }
# }


compute_res(trainingDays, 50, FALSE, training_days_method)
