# code to compute test statistics of DM tests
# results are saved in a specific data frame format to simplify plotting later on

# rm(list=ls())


fix_training_days <- FALSE
training_days_method <- "last_m_days"

compute_dm <-function(timeWindow, fix_training_days, training_days_method) {

  for (groupNR in 1:5) {
    benchmark <- "ssh.i"
    
    
    fName <- paste0("Res_group_", groupNR)
    
    if (fix_training_days) {
      fName <- paste0(training_days_method, "_", fName)
    }
    
    fName <- paste0("m_", timeWindow, "_", fName)
    
    load(paste0("Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable
    
    print(fName)
    
    
    library(forecast) # for DM test function
    
    input_models <- names(res$es_list)
    input_models <- input_models[! input_models %in% c(benchmark, "ens","emos.q","ecc.q","ecc.s","decc.q")]
    input_scores <- names(res)
    input_scores <- input_scores[!input_scores %in% c("param_list", "indep_list", "tau", "timing_list", "chosenCopula_list", "obs", "mvpp_list")]
    input_scores <- c(input_scores, "crps_1", "crps_2", "crps_3")
    
    dfmc <- expand.grid(input_models, input_scores, NA)
    names(dfmc) <- c("model", "score", "value")
      
    options(warn=1)
    
    for(this_model in input_models){
      for(this_score in input_scores){
        ind <- which(dfmc$model == this_model & dfmc$score == this_score)  
        
        # deal with CRPS specifically (use only first dimension, not all 5 recorded ones)
        
        # print(paste0("( ", this_model, " , ", this_score, ")"))
        if(this_score %in% c("crps_1", "crps_2", "crps_3")){
          
          n <- as.numeric(gsub("[^0-9]", "", this_score))
          
          tmp <- NA
          tryDM <- try(tmp_DM <- dm.test(e1 = res[["crps_list"]][[this_model]][, n],
                                         e2 = res[["crps_list"]][[benchmark]][, n],
                                         h = 1, power = 1), silent = FALSE)
          if(class(tryDM) != "try-error"){
            tmp <- tmp_DM$statistic
          } else{
            tmp <- 0
          }
          dm_teststat_vec <- tmp
        }else if(this_score == "crps_list"){
    
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
        
        # print(dm_teststat_vec)
        # print("")
        
        
        dfmc$value[ind] <- dm_teststat_vec
      }
    }
    
    savedir <- paste0("Data/TestStatistic/")
    savename <- paste0("TestStatistic_data_group",groupNR, ".Rdata")
    
    if (fix_training_days) {
      savename <- paste0(training_days_method, "_", savename)
    }
    

    savename <- paste0("m_", timeWindow, "_", savename)

  
    save(dfmc, file = paste0(savedir, savename))
  }
}

# for (fix_training_days in c(TRUE, FALSE)) {
#   if (fix_training_days) {
#     for (training_days_method in c("random_past", "last_m_days", "random_2w_interval")) {
#       
#       compute_dm(100, fix_training_days, training_days_method)
#       
#     }
#   } else {
#     compute_dm(100, fix_training_days, training_days_method)
#   }
# }

# Standard time window of 50
compute_dm(50, FALSE, training_days_method)
