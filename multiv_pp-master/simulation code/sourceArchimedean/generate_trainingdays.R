getCentralDays <- function(value, minRange, maxRange) {
  DAYS_PER_YEAR <- 365
  
  # Find smallest value
  smallest_value <- as.integer(value - DAYS_PER_YEAR * floor((value - minRange) / DAYS_PER_YEAR))
  
  # Find largest value
  largest_value <- as.integer(value + DAYS_PER_YEAR * floor((maxRange - value) / DAYS_PER_YEAR))
  
  return(seq(smallest_value, largest_value, by = DAYS_PER_YEAR))
}

getIntervals <- function(day, minRange, maxRange, interval_length) {
  central_days <- getCentralDays(day, minRange, maxRange)
  
  intervals <- unlist(lapply(central_days, function(x) {
    seq(max(minRange, x - interval_length), min(maxRange, x + interval_length))
  }))
  
  return(intervals)
}

getTrainingDays <- function(initial_days, n, m, method) {
  
  set.seed(1)
  
  out <- data.frame(TrainingDays = I(vector("list", n)))
  
  if (method == "random_past") {
    for (nn in 1:n) {
      
      obs_IDs <- sample(x = 1:(initial_days + nn - 1), size = m, replace = FALSE)
      
      out$TrainingDays[[nn]] <- obs_IDs
      
    }
  } else if (method == "last_m_days") {
    for (nn in 1:n) {
      
      obs_IDs <- (initial_days + nn - m):(initial_days + nn - 1)
      
      out$TrainingDays[[nn]] <- obs_IDs
      
    }
  } else if (method == "random_2w_interval") {
    for (nn in 1:n) {
      
      obs_IDs <- sample(x = getIntervals(nn, 1, n, 28), size = m, replace = FALSE)
      
      out$TrainingDays[[nn]] <- obs_IDs
      
    }
  }
  
  return(out)
}