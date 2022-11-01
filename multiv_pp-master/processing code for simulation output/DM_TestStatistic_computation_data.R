# code to compute test statistics of DM tests
# results are saved in a specific data frame format to simplify plotting later on

# rm(list=ls())

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/University/Bachelor/Year 3/BEP - BAM/Code/multiv_pp-master/processing code for simulation output")



groupNR <- 3
fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable



library(forecast) # for DM test function

input_models <- names(res$es_list)
input_scores <- names(res)

dfmc <- expand.grid(input_models, input_scores[!input_scores %in% c("param_list", "indep_list", "tau", "timing_list")], NA)
names(dfmc) <- c("model", "score", "value")
  
for(this_model in input_models){
  for(this_score in input_scores){
    ind <- which(dfmc$model == this_model & dfmc$score == this_score)  
    
    # deal with CRPS specifically (use only first dimension, not all 5 recorded ones)
    
    if(this_score == "crps_list"){

      tmp <- NA
      tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]],
                                     e2 = res[[which(input_scores == this_score)]][[which(input_models == "gca")]],
                                     h = 1, power = 1), silent = TRUE)
      if(class(tryDM) != "try-error"){
        tmp <- tmp_DM$statistic
      } else{
        tmp <- 0
      }
      dm_teststat_vec <- tmp
    } else{
        tmp <- NA
        tryDM <- try(tmp_DM <- dm.test(e1 = res[[which(input_scores == this_score)]][[which(input_models == this_model)]],
                                       e2 = res[[which(input_scores == this_score)]][[which(input_models == "gca")]],
                                       h = 1, power = 1), silent = TRUE)
        if(class(tryDM) != "try-error"){
          tmp <- tmp_DM$statistic
        } else{
          tmp <- 0
        }
        dm_teststat_vec <- tmp
    }
    
    
    dfmc$value[ind] <- dm_teststat_vec
  }
}



save(dfmc, file = paste0("../Data/TestStatistic/TestStatistic_data_group",groupNR, ".Rdata"))

