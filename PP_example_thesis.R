
# round(rnorm(5, mean = 285, sd = 1),2)


# Observations
obs_one <- c(296.80, 292.85, 299.05, 287.76, 289.68)
obs_two <- c(284.90, 284.97, 285.06, 285.58, 285.29)

# ensemble members
ens_one <- c(284.70, 282.20, 292.68, 295.75, 284.78)
ens_two <- c(284.50, 283.93, 285.74, 282.70, 284.66)

7+mean(obs_one)

moments <- function(x, m) {
  return(sum(x^m)/length(x))
}

moments(obs_one, 1)
moments(obs_one, 1)


moments(obs_one, 2)
moments(obs_one, 2)


sdestimate <- function(x) {
  return(sqrt(moments(x, 2) - moments(x,1)^2))
}

sdestimate(obs_one)^2
sdestimate(obs_one)^2

toTable <- function(x,y) {
  textresult <- ""
  for (i in 1:length(ens_one)) {
    textresult <- paste0(textresult, i, " & (", rank(x)[i], ") ", x[i], " & (", rank(y)[i], ") ", y[i], "\\ ")
  }
  return(textresult)
}

set.seed(1)
sample_one <- round(rnorm(5, mean = 300.23, sd = sqrt(20.12)),2)
sample_two <- round(rnorm(5, mean = 285.16, sd = sqrt(0.06)),2)

toTable(sample_one, sample_two)

indexShift <- function(x,y) {
  index <- c()
  for (i in 1:5) {
    index <- c(index, match(rank(x)[i], rank(y))[1])
  }
  return(index)
}

shuffle_sample <- function(x,y) {
  index <- indexShift(x,y)
  out <- array(NA, dim=5)
  for (i in 1:5) {
    out[index[i]] <- x[i]
  }
  return(out)
}
indexShift(sample_one, ens_one)
shuffle_sample(sample_one, ens_one)
toTable(shuffle_sample(sample_one, obs_one), shuffle_sample(sample_two, obs_two))
