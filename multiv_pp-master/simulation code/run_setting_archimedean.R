rm(list=ls())


library(scoringRules)
library(MASS)

# "source" directory 
setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/simulation code")
dir <- "./sourceArchimedean/"
source(paste0(dir, "generate_observations_arch.R"))
source(paste0(dir, "generate_ensfc_arch.R"))
source(paste0(dir, "postprocess_ensfc_arch.R"))
source(paste0(dir, "mvpp_arch.R"))
source(paste0(dir, "evaluation_functions_arch.R"))


# Any change in the number of parameters should be reflected in the run_wrapper function

# Dependence parameter between weather variables for observations
input_theta0 <- c(5, 10)
# Dependence parameter between weather variables for ensemble forecasts
input_theta <- c(5, 10)

# Copula type
input_copula <- c("Frank","Gumbel", "Clayton")

# Dimension of multi-index l (weather variable, position and look-ahead time)
input_d <- 3

# Put the parameters in a grid
input_par <- expand.grid(input_theta0, input_theta, input_copula, input_d)
names(input_par) <- c("theta0", "theta", "copula",  "d")

# Number of Monte Carlo repetitions
MC_reps <- 10

# Parameter that influences reliability of methods such as Ssh and ECC.S
randomRepetitions <- 5


# Save data
Rdata_dir <- "../Data/Rdata/Archimedean" # directory to save Rdata files to
Rout_dir <- "../Data/Rout/Archimedean" # directory to save Rout files to


# Simulation parameters
evalDays <- 2
trainingDays <- 5
ensembleMembers <- 3

# Model selection. 
# ObservationsModel denotes the model that is used to generate the observations
# ForecastModel denotes the model that is used to generate the forecasts
observationsModel <- 1
forecastModel <- 1




eval_all_mult <- function(mvpp_out, obs){
  esout <- es_wrapper(mvpp_out, obs)
  vs1out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 1)
  vs1wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 1)
  vs0out <- vs_wrapper(mvpp_out, obs, weight = FALSE, p = 0.5)
  vs0wout <- vs_wrapper(mvpp_out, obs, weight = TRUE, p = 0.5)
  return(list("es" = esout, "vs1" = vs1out, "vs1w" = vs1wout, "vs0" = vs0out, "vs0w" = vs0wout))
}

run_setting1 <- function(obsmodel, fcmodel, nout, ninit, nmembers, d, MCrep, rand_rep, progress_ind = FALSE, compute_crps, ...){
  
  # generate objects to save scores to
  modelnames <- c("ens", "emos.q", "ecc.q", "ecc.s", "decc.q", "ssh", "gca","clayton","frank","gumbel")
  crps_list <- es_list <- vs1_list <- vs1w_list <- vs0_list <- vs0w_list <- list()
  for(mm in 1:length(modelnames)){
    es_list[[mm]] <- vs1_list[[mm]] <- vs1w_list[[mm]] <- 
      vs0_list[[mm]] <- vs0w_list[[mm]] <- matrix(NA, nrow = MCrep, ncol = nout) 
    crps_list[[mm]] <- array(NA, dim = c(MCrep, nout, d))
  }
  names(crps_list) <- names(es_list) <- names(vs1_list) <- names(vs1w_list) <- 
    names(vs0_list) <- names(vs0w_list) <- modelnames
  
  for(rr in 1:MCrep){
    if(progress_ind){
      if(rr %% 1 == 0){
        cat("starting at", paste(Sys.time()), ": MC repetition", rr, "of", MCrep, "\n"); flush(stdout())
      }
    }
    
    # set random seed
    set.seed(42+rr)
    
    # generate observations
    obs <- generate_obs(model = obsmodel, nout = nout, ninit = ninit, d = d, ...)
    
    # generate ensemble forecasts
    fc <- generate_ensfc(model = fcmodel, nout = nout, ninit = ninit, nmembers = nmembers, d = d, ...)
    
    # postprocess ensemble forecasts
    pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init, 
                       obs = obs$obs, obs_init = obs$obs_init, 
                       train = "init", trainlength = NULL, emos_plus = TRUE)
    
    # if there are NaN's in pp output, generate new sets of obs and fc, and re-run pp code
    # only happens very rarely for extreme parameter combinations
    while(anyNA(pp_out)){
      obs <- generate_obs(model = obsmodel, nout = nout, ninit = ninit, d = d, ...)
      fc <- generate_ensfc(model = fcmodel, nout = nout, ninit = ninit, nmembers = nmembers, d = d, ...)
      pp_out <- postproc(fcmodel = fcmodel, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init, 
                         obs = obs$obs, obs_init = obs$obs_init, 
                         train = "init", trainlength = NULL, emos_plus = TRUE)
    }
    
    ## iterate over models and compute scores

    # ensemble forecasts
    if(compute_crps){
      crps_list$ens[rr, , ] <- crps_wrapper(fc$ensfc, obs$obs)
    }
    tmp <- eval_all_mult(mvpp_out = fc$ensfc, obs = obs$obs)
    es_list$ens[rr, ] <- tmp$es 
    vs1_list$ens[rr, ] <- tmp$vs1
    vs1w_list$ens[rr, ] <- tmp$vs1w
    vs0_list$ens[rr, ] <- tmp$vs0
    vs0w_list$ens[rr, ] <- tmp$vs0w
    
    # EMOS.Q
    emos.q <- mvpp(method = "EMOS", variant = "Q", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                   obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
    
    if(compute_crps){
      crps_list$emos.q[rr, , ] <- crps_wrapper(emos.q, obs$obs)
    }
    
    tmp <- eval_all_mult(mvpp_out = emos.q, obs = obs$obs)
    es_list$emos.q[rr, ] <- tmp$es 
    vs1_list$emos.q[rr, ] <- tmp$vs1
    vs1w_list$emos.q[rr, ] <- tmp$vs1w
    vs0_list$emos.q[rr, ] <- tmp$vs0
    vs0w_list$emos.q[rr, ] <- tmp$vs0w
    
    # ECC.Q
    ecc.q <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, 
                  EMOS_sample = emos.q)
    
    if(compute_crps){
      crps_list$ecc.q[rr, , ] <- crps_wrapper(ecc.q, obs$obs)
    }
    tmp <- eval_all_mult(mvpp_out = ecc.q, obs = obs$obs)
    es_list$ecc.q[rr, ] <- tmp$es 
    vs1_list$ecc.q[rr, ] <- tmp$vs1
    vs1w_list$ecc.q[rr, ] <- tmp$vs1w
    vs0_list$ecc.q[rr, ] <- tmp$vs0
    vs0w_list$ecc.q[rr, ] <- tmp$vs0w
    
    # ECC.S -> involves randomness -> repeat rand_rep times
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    for(RR in 1:rand_rep){
      emos.s <- mvpp(method = "EMOS", variant = "S", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                     obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out)
      ecc.s <- mvpp(method = "ECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, 
                    EMOS_sample = emos.s)
      
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(ecc.s, obs$obs)
      }
      tmp <- eval_all_mult(mvpp_out = ecc.s, obs = obs$obs)
      es_list_tmp[,RR] <- tmp$es
      vs1_list_tmp[,RR] <- tmp$vs1
      vs1w_list_tmp[,RR] <- tmp$vs1w
      vs0_list_tmp[,RR] <- tmp$vs0
      vs0w_list_tmp[,RR] <- tmp$vs0w
    }
    crps_list$ecc.s[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
    es_list$ecc.s[rr, ] <- apply(es_list_tmp, 1, mean)
    vs1_list$ecc.s[rr, ] <- apply(vs1_list_tmp, 1, mean)
    vs1w_list$ecc.s[rr, ] <- apply(vs1w_list_tmp, 1, mean)
    vs0_list$ecc.s[rr, ] <- apply(vs0_list_tmp, 1, mean)
    vs0w_list$ecc.s[rr, ] <- apply(vs0w_list_tmp, 1, mean)
    
    # dECC.Q
    decc.q <- mvpp(method = "dECC", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                   obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out, 
                   EMOS_sample = emos.q, ECC_out = ecc.q) 
    
    if(compute_crps){
      crps_list$decc.q[rr, , ] <- crps_wrapper(decc.q, obs$obs)
    }
    tmp <- eval_all_mult(mvpp_out = decc.q, obs = obs$obs)
    es_list$decc.q[rr, ] <- tmp$es 
    vs1_list$decc.q[rr, ] <- tmp$vs1
    vs1w_list$decc.q[rr, ] <- tmp$vs1w
    vs0_list$decc.q[rr, ] <- tmp$vs0
    vs0w_list$decc.q[rr, ] <- tmp$vs0w
    
    # SSh -> involves randomness -> repeat rand_rep times
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    for(RR in 1:rand_rep){
      ssh <- mvpp(method = "SSh", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out,
                  EMOS_sample = emos.q) 
      
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(ssh, obs$obs)
      }
      tmp <- eval_all_mult(mvpp_out = ssh, obs = obs$obs)
      es_list_tmp[,RR] <- tmp$es
      vs1_list_tmp[,RR] <- tmp$vs1
      vs1w_list_tmp[,RR] <- tmp$vs1w
      vs0_list_tmp[,RR] <- tmp$vs0
      vs0w_list_tmp[,RR] <- tmp$vs0w
    }
    crps_list$ssh[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
    es_list$ssh[rr, ] <- apply(es_list_tmp, 1, mean)
    vs1_list$ssh[rr, ] <- apply(vs1_list_tmp, 1, mean)
    vs1w_list$ssh[rr, ] <- apply(vs1w_list_tmp, 1, mean)
    vs0_list$ssh[rr, ] <- apply(vs0_list_tmp, 1, mean)
    vs0w_list$ssh[rr, ] <- apply(vs0w_list_tmp, 1, mean)
    
    # GCA -> involves randomness -> repeat rand_rep times
    es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
      vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
    crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
    for(RR in 1:rand_rep){
      gca <- mvpp(method = "GCA", ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                  obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out) 
      
      if(compute_crps){
        crps_list_tmp[,,RR] <- crps_wrapper(gca, obs$obs)
      }
      tmp <- eval_all_mult(mvpp_out = gca, obs = obs$obs)
      es_list_tmp[,RR] <- tmp$es
      vs1_list_tmp[,RR] <- tmp$vs1
      vs1w_list_tmp[,RR] <- tmp$vs1w
      vs0_list_tmp[,RR] <- tmp$vs0
      vs0w_list_tmp[,RR] <- tmp$vs0w
    }
    crps_list$gca[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
    es_list$gca[rr, ] <- apply(es_list_tmp, 1, mean)
    vs1_list$gca[rr, ] <- apply(vs1_list_tmp, 1, mean)
    vs1w_list$gca[rr, ] <- apply(vs1w_list_tmp, 1, mean)
    vs0_list$gca[rr, ] <- apply(vs0_list_tmp, 1, mean)
    vs0w_list$gca[rr, ] <- apply(vs0w_list_tmp, 1, mean)
    
    
    # Archimedean copulas -> involves randomness -> repeat rand_rep times
    for (method in c("Clayton","Frank", "Gumbel")) {
      es_list_tmp <- vs1_list_tmp <- vs1w_list_tmp <- 
        vs0_list_tmp <- vs0w_list_tmp <- matrix(NA, nrow = nout, ncol = rand_rep)
      crps_list_tmp <- array(NA, dim = c(nout, d, rand_rep))
      for(RR in 1:rand_rep){
        mvd <- mvpp(method = method, ensfc = fc$ensfc, ensfc_init = fc$ensfc_init,
                    obs = obs$obs, obs_init = obs$obs_init, postproc_out = pp_out) 
        
        if(compute_crps){
          crps_list_tmp[,,RR] <- crps_wrapper(mvd, obs$obs)
        }
        tmp <- eval_all_mult(mvpp_out = mvd, obs = obs$obs)
        es_list_tmp[,RR] <- tmp$es
        vs1_list_tmp[,RR] <- tmp$vs1
        vs1w_list_tmp[,RR] <- tmp$vs1w
        vs0_list_tmp[,RR] <- tmp$vs0
        vs0w_list_tmp[,RR] <- tmp$vs0w
      }
      
      if (method == "Clayton") {
        crps_list$clayton[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
        es_list$clayton[rr, ] <- apply(es_list_tmp, 1, mean)
        vs1_list$clayton[rr, ] <- apply(vs1_list_tmp, 1, mean)
        vs1w_list$clayton[rr, ] <- apply(vs1w_list_tmp, 1, mean)
        vs0_list$clayton[rr, ] <- apply(vs0_list_tmp, 1, mean)
        vs0w_list$clayton[rr, ] <- apply(vs0w_list_tmp, 1, mean)
        print("DONE")
      } else if (method == "Frank") {
        crps_list$frank[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
        es_list$frank[rr, ] <- apply(es_list_tmp, 1, mean)
        vs1_list$frank[rr, ] <- apply(vs1_list_tmp, 1, mean)
        vs1w_list$frank[rr, ] <- apply(vs1w_list_tmp, 1, mean)
        vs0_list$frank[rr, ] <- apply(vs0_list_tmp, 1, mean)
        vs0w_list$frank[rr, ] <- apply(vs0w_list_tmp, 1, mean)
      } else if (method == "Gumbel") {
        crps_list$gumbel[rr,,] <- apply(crps_list_tmp, c(1,2), mean)
        es_list$gumbel[rr, ] <- apply(es_list_tmp, 1, mean)
        vs1_list$gumbel[rr, ] <- apply(vs1_list_tmp, 1, mean)
        vs1w_list$gumbel[rr, ] <- apply(vs1w_list_tmp, 1, mean)
        vs0_list$gumbel[rr, ] <- apply(vs0_list_tmp, 1, mean)
        vs0w_list$gumbel[rr, ] <- apply(vs0w_list_tmp, 1, mean)
      } 
    }
    
    # end loop over Monte Carlo repetitions  
  }
  
  # return results, as a huge list
  out <- list("crps_list" = crps_list, "es_list" = es_list, "vs1_list" = vs1_list,
              "vs1w_list" = vs1w_list, "vs0_list" = vs0_list, "vs0w_list" = vs0w_list)
  return(out)
}


run_wrapper <- function(runID){
  sink(file = paste0(Rout_dir, "setting1_", runID, ".Rout"))
  par_values <- as.numeric(input_par[ID, ])
  res <- run_setting1(obsmodel = observationsModel, fcmodel = forecastModel, nout = evalDays, ninit = trainingDays, nmembers = ensembleMembers, 
                  MCrep = MC_reps, rand_rep = randomRepetitions, 
                  progress_ind = TRUE, compute_crps = TRUE,
                  theta0 = input_par$theta[runID], 
                  theta = input_par$theta[runID], 
                  copula = input_par$copula[runID],
                  d = input_par$d[runID])
  savename <- paste0(Rdata_dir, "res_setting_arch", runID, ".Rdata")
  save(res, input_par, file = savename)
  sink()
}


# Run the simulation for all parameter coefficients with a unique ID
for (ID in 1:dim(input_par)[1]) {
  closeAllConnections()
  print(ID)
  run_wrapper(runID = ID)
}


