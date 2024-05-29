# code to compute test statistics of DM tests
# results are saved in a specific data frame format to simplify plotting later on

# rm(list=ls())

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

for (groupNR in 1:5) {
  benchmark <- "ssh.h"
  
  fix_training_days <- FALSE
  training_days_method <- "last_m_days"
  
  
  fName <- paste0("Res_group_", groupNR)
  
  if (fix_training_days) {
    fName <- paste0(training_days_method, fName)
  }
  
  load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable
  
  
  
  library(forecast) # for DM test function
  
  input_models <- names(res$es_list)
  input_models <- input_models[! input_models %in% c(benchmark)]
  input_scores <- names(res)
  input_scores <- input_scores[!input_scores %in% c("param_list", "indep_list", "tau", "timing_list", "chosenCopula_list", "obs", "mvpp_list")]
  
  dfmc <- expand.grid(input_models, input_scores, NA)
  names(dfmc) <- c("model", "score", "value")
    
  options(warn=1)
  
  for(this_model in input_models){
    for(this_score in input_scores){
      ind <- which(dfmc$model == this_model & dfmc$score == this_score)  
      
      # deal with CRPS specifically (use only first dimension, not all 5 recorded ones)
      
      print(paste0("( ", this_model, " , ", this_score, ")"))
      
      if(this_score == "crps_list"){
  
        tmp <- NA
        tryDM <- try(tmp_DM <- dm.test(e1 = res[[this_score]][[this_model]],
                                       e2 = res[[this_score]][[benchmark]],
                                       h = 1, power = 1), silent = FALSE)
        if(class(tryDM) != "try-error"){
          tmp <- tmp_DM$statistic
        } else{
          tmp <- 0
        }
        dm_teststat_vec <- tmp
      } else{
          tmp <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[[this_score]][[this_model]],
                                         e2 = res[[this_score]][[benchmark]],
                                         h = 1, power = 1), silent = FALSE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
          } else{
            tmp <- 0
          }
          dm_teststat_vec <- tmp
      }
      
      print(dm_teststat_vec)
      print("")
      
      
      dfmc$value[ind] <- dm_teststat_vec
    }
  }
  
  savedir <- paste0("../Data/TestStatistic/")
  savename <- paste0("TestStatistic_data_group",groupNR, ".Rdata")
  
  if (fix_training_days) {
    savename <- paste0(training_days_method, savename)
}

save(dfmc, file = paste0(savedir, savename))
}
