library(MASS)
library(copula)

set.seed(2)

# Generate data
n <- 50
S <- matrix(c(2, 0.5, 0.2,
              0.5, 1, 0.3,
              0.2, 0.3, 3), 
            nrow = 3, ncol = 3, byrow = TRUE)
mu <- c(X = 50, Y = -3, Z = 7)
sigma <- sqrt(diag(S))
data <- mvrnorm(n, mu = mu, Sigma = S)

d <- dim(S)[1]

obs_train_CDF <- c()
for (dd in 1:d) {
  obs_train_CDF <- cbind(obs_train_CDF, pnorm(data[,dd], mean = mu[dd], sd = sigma[dd]))
}

print("-- Original data --")
print("Covariance matrix:")
print(S)
print("SSE data vs ground truth:")
print(sum((cov(data) - S)^2))
print("Mean vector:")
print(mu)
print("")

####################
## COPULA Package ##
####################

# Fit normal copula
cop <- fitCopula(normalCopula(dim = 3, dispstr = "un"), data = obs_train_CDF, method = "mpl", optim.control = list(maxit = 1000))

params <- mapply(function(mean, sd) list(mean = mean, sd = sd), mu, sigma, SIMPLIFY = FALSE)

mvDistribution <- mvdc(copula=cop@copula, margins=rep("norm", d),
                       paramMargins=params)

X_copula <- rMvdc(n, mvDistribution)

print("-- Copula package --")
print("Covariance matrix:")
print(cov(X_copula))
print("Correlation matrix:")
print(cor(X_copula))
print("SSE:")
print(sum((cov(X_copula) - S)^2))
print("Mean vector:")
print(apply(X_copula, 2, mean))
print("")

#####################
## Manual approach ##
#####################

obs_values <- c()
for (dd in 1:d) {
  obs_values <- cbind(obs_values, qnorm(obs_train_CDF[,dd]))
}

cov_obs <- cov(obs_values)

# Dependence structure by Copula
Z_manual <- mvrnorm(n = n, mu = rep(0, d), Sigma = cov_obs)

X_manual <- c()
for (dd in 1:d) {
  X_manual <- cbind(X_manual, Z_manual[,dd] * sigma[dd] + mu[dd])
}

print("-- Manual approach --")
print("Covariance matrix:")
print(cov(X_manual))
print("SSE:")
print(sum((cov(X_manual) - S)^2))
print("Mean vector:")
print(apply(X_manual, 2, mean))
print("")

