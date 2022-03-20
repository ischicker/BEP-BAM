library(MASS)
library(psych)

# Reproducibility
set.seed(100)

# Numer of variables
m <- 3

# Number of observations
n <- 2000

# Correlation matrix
sigma <- matrix(c(1, 0.4, 0.2,
                  0.4, 1, -0.8,
                  0.2, -0.8, 1), 
                nrow=3)

# Create multivariate observations
z <- mvrnorm(n,mu=rep(0, m),Sigma=sigma,empirical=T)

# Estimate the correlation matrix
cor(z,method='spearman')

# Create scatter plots, linear regressions and visualisations of the data
pairs.panels(z)

# Since z is distributed with as a standard normal, pnorm(z) yields a uniform distribution 
u <- pnorm(z)
pairs.panels(u)

# Now using the copula package
library(copula)
set.seed(100)

# Define the dependence structure with the copula
myCop <- normalCopula(param=c(0.4,0.2,-0.8), dim = 3, dispstr = "un")

# Construct the multivariate distribution (Constructed from copula)
myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta", "t"),
              paramMargins=list(list(shape=2, scale=1),
                                list(shape1=2, shape2=2), 
                                list(df=5)) )

# Random sample from the mvdc
Z2 <- rMvdc(2000, myMvd)
colnames(Z2) <- c("x1", "x2", "x3")
pairs.panels(Z2)
