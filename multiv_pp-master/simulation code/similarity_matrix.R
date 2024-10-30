
source("simulation code/similarity_criterion.R")

create_ensfc <- function(data) {
  
  trainingDays <- 365
  
  stations <- unique(data$stat)
  days <- sort(unique(data$td))
  nout <- length(days) - trainingDays
  
  # 1 weather variable over 1 lead time with d stations
  d <- length(stations)
  
  # Ensemble members
  m <- sum(grepl("laef", names(data)))
  ensembleMembers <- sapply(1:m, FUN = function(x) paste0("laef", x))
  
  # Create ensfc and obs data structure
  ensfc <- array(NA, dim = c(nout, m, d))

  
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
          if (day > trainingDays) {
            ensfc[day - trainingDays, , index] <- unlist(dat[ensembleMembers], use.names = FALSE)
          }
          
        }
        ,
        error = function(cond) {
          print(cond)
        }
      )
    }
  }
  
  return (ensfc)
}

compute_sim_matrix <- function(groupNr) {
  
  # Load the input data
  x <- load(paste0("Data/Input/data", groupNr, ".Rdata"))
  assign("data", get(x))
  rm(list = x, envir = .GlobalEnv)
  
  x <- load(paste0("Data/UVPP/uvpp", groupNr, ".Rdata"))
  assign("uvpp", get(x))
  rm(list = x, envir = .GlobalEnv)
  
  ensfc <- create_ensfc(data)

  # Vectorize sim
  sim_vectorized <- Vectorize(function(i, j) Sim(i, j, uvpp, ensfc))
  
  # Generate indices for upper triangular part
  n <- dim(ensfc)[1]
  indices <- combn(n, 2)
  
  # Compute similarity scores for all pairs
  sim_values <- sim_vectorized(indices[1, ], indices[2, ])
  
  # Initialize the similarity matrix
  sim_matrix <- matrix(0, n, n)
  
  # Fill the similarity matrix
  sim_matrix[upper.tri(sim_matrix)] <- sim_values
  sim_matrix <- sim_matrix + t(sim_matrix)
  
  # Save results
  savename <- paste0("Data/SimilarityMatrix/simMatrix", groupNr, ".Rdata")
  save(sim_matrix, file = savename)

}

# for (groupNr in 1:5) {
#   compute_sim_matrix(groupNr)
# }
