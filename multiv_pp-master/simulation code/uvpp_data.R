
compute_uvpp <- function(groupNr) {
  
  # Load the input data
  x <- load(paste0("Data/Input/data", groupNr, ".Rdata"))
  assign("data", get(x))
  rm(list = x, envir = .GlobalEnv)
  
  uvpp <- list()
  stations <- unique(data$stat)
  trainingDays <- 365
  
  for (station in stations) {
    
    statIndex <- match(station, stations)
    
    emos.result <- emos_T2M_mean_singleForecast(subset(data, stat == station), trainingDays)
    emos.result$stat <- statIndex
    
    uvpp <- rbind(uvpp, emos.result)
    
  }

  savename <- paste0("Data/UVPP/uvpp", groupNr, ".Rdata")
  save(uvpp, file = savename)
}

# for (groupNr in 1:5) {
#   compute_uvpp(groupNr)
# }