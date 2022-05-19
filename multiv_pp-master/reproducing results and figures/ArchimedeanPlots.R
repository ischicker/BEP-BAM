require(MASS)
require(copula)
require(ggplot2)
require(FactorCopula)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")


tau <- 0.8
d <- 2

dfplot <- data.frame(matrix(ncol=3, nrow=0))
names(dfplot) <- c("copula", "xvalue", "yvalue")
cops <- c("Frank", "Gumbel", "Clayton")

n <- 5000



for (copula in cops) {
  # Generate observations
  if (copula == "Frank"){
    # -infty < theta0 < infty and theta0 != 0
    param <- tau2par("frk", tau)
    cop <- frankCopula(param = param, dim = d)
  } else if (copula == "Gumbel") {
    # 1 <= theta0 < infty
    param <- tau2par("gum", tau)
    cop <- gumbelCopula(param = param, dim = d)
  } else if (copula == "Clayton") {
    # -1 <= theta0 < infty and theta0 != 0
    param <- 2 * tau / (1-tau)
    cop <- claytonCopula(param = param, dim = d)
  } else {
    stop(paste("Incorrect value for 'copula' input"))
  }
  
  paramMargins <- list()
  
  for (i in 1:d){
    paramMargins[[i]] <- list(mean = 0, sd = 1)
  }
  
  # Generate observations with standard normal marginals
  mvDistribution <- mvdc(copula=cop, margins=rep("norm", d),
                         paramMargins=paramMargins)
  
  # Random sample from the mvdc
  obs <- rMvdc(n, mvDistribution)
  
  new_df <- data.frame("copula" = copula, "xvalue" = obs[,1], "yvalue" = obs[,2])
  dfplot <- rbind(dfplot, new_df)

}

p1 <- ggplot(dfplot, aes(x=xvalue, y=yvalue)) + geom_point()

p1 <- p1 + facet_grid(cols = vars(copula),
                      labeller = label_bquote(cols = copula: .(as.character(copula))))


p1 <- p1 + theme_bw()
p1 <- p1 + xlab(bquote(x)) + ylab(bquote(y)) + coord_fixed(ratio=1)


fileName <- paste0("../Data/Plots/ArchimedeanPlots.png")

plotWidth <- 9
plotHeight <- 6
res <- 400


ggsave(
  fileName,
  p1,
  width = plotWidth,
  height = plotHeight,
  dpi = res
)

