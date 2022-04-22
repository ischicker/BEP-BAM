# code to compute test statistics of DM tests
# results are saved in a specific data frame format to simplify plotting later on

rm(list=ls())

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/processing code for simulation output")

# parameters 
# Should be the same as in run_setting?.R
# parameters to run

# Dependence parameter between weather variables for observations
# input_theta0 <- c(5, 10)
# Dependence parameter between weather variables for ensemble forecasts
# input_theta <- c(5, 10)

# Copula type
input_copula <- c("Frank","Gumbel", "Clayton")

# Dimension
input_d <- 3

# Repetitions
repetitions <- 10

# Put the parameters in a grid
input_par <- expand.grid(input_copula, input_d, 1:repetitions)
names(input_par) <- c("copula",  "d", "repetition")

# Number of Monte Carlo repetitions
MC_reps <- 75

# Setting parameter for different runs
setting <- 1

# Model 1 : Standard Gaussian Marginals

observationsModel <- 2


forecastModel <- 2


fName <- paste0("Archimedean","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,"_ID_")

df_raw <- data.frame(input_par)
df_raw$simID <- 1:nrow(df_raw)

flist <- list.files("../Data/Rdata/")
existing <- as.numeric(sapply(flist, FUN = function(x) as.numeric(strsplit(strsplit(x, fName)[[1]][2], ".Rdata"))))



df <- df_raw[which(is.element(df_raw$simID, existing)),] 

load(paste0("../Data/Rdata/",fName, "1.Rdata"))

input_models <- names(res$es_list)
input_scores <- names(res)

df_use <- as.data.frame(df[1,])
df_use$model <- as.character("a")
df_use$score <- as.character("a")

model_score_grid <- expand.grid(input_models, input_scores[!input_scores %in% c("param_list", "indep_list", "tau", "timing_list")])

# Produces df of the form
#         rho0  eps   sigma   rho   d   simID   model   score
#   1     0.1   0     0.5     0.1   3   1       ens     crps_list
for(i in 1:nrow(df)){
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),] <- df[i,]
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$model <- as.character(model_score_grid$Var1)
  df_use[((i-1)*nrow(model_score_grid)+1):(i*nrow(model_score_grid)),]$score <- as.character(model_score_grid$Var2)
}

# Makes a copy of df_use for each MC_rep
dfmc <- data.frame(cbind(zoo::coredata(df_use)[rep(seq(nrow(df_use)),MC_reps),]))
dfmc$value <- NA
dfmc$tau <- NA
# head(dfmc)

library(forecast) # for DM test function

for(ID in existing[!is.na(existing)]){
  load(paste0("../Data/Rdata/", fName, ID,".Rdata"))
  print(ID)

  
  for(this_model in input_models){
    for(this_score in input_scores){
      ind <- which(dfmc$simID == ID & dfmc$model == this_model & dfmc$score == this_score)  
      
      # deal with CRPS specifically (use only first dimension, not all 5 recorded ones)
      
      if(this_score == "crps_list"){
        dm_teststat_vec <- rep(NA, MC_reps)
        for(MC_rep in 1:MC_reps){
          tmp <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]][,,1][MC_rep,],
                                         e2 = res[[which(input_scores == this_score)]][[which(input_models == "ecc.q")]][,,1][MC_rep,],
                                         h = 1, power = 1), silent = TRUE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
          } else{
            tmp <- 0
          }
          dm_teststat_vec[MC_rep] <- tmp
        }
      } else{
        dm_teststat_vec <- rep(NA, MC_reps)
        for(MC_rep in 1:MC_reps){
          tmp <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]][MC_rep,],
                                         e2 = res[[which(input_scores == this_score)]][[which(input_models == "ecc.q")]][MC_rep,],
                                         h = 1, power = 1), silent = TRUE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
          } else{
            tmp <- 0
          }
          dm_teststat_vec[MC_rep] <- tmp
        }
      }
      
      
      dfmc$value[ind] <- dm_teststat_vec
      dfmc$tau[ind] <- res$tau
    }
  }
  
}


save(dfmc, file = paste0("../Data/TestStatistic/TestStatistic_Archimedean","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".Rdata"))

