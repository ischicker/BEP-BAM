groupNR <- 1
fix_training_days <- TRUE
training_days_method <- "random_2w_interval" #"random_past" # "random_2w_interval" # "last_m_days"
timeWindow <- 0
    
fName <- paste0("Res_group_", groupNR)

if (fix_training_days) {
  if (timeWindow == 100){
  fName <- paste0(training_days_method, "_", fName)
  } else {
    fName <- paste0(training_days_method, fName)
  }
}

if (timeWindow == 100){
  fName <- paste0("m_", timeWindow, "_", fName)
}

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

load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable

differences <- sum(res$mvpp_list$ssh.h - res$mvpp_list$ssh.i != 0)
total_values <- prod(dim(res$mvpp_list$ssh.h))
percentage <- differences / total_values * 100
print("Percentage of collisions:")
print(percentage)
