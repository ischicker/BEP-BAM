require(MASS)
require(copula)
require(FactorCopula)
require(ggplot2)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")


MC_reps <- 100
tau <- 0.5
d <- 3
days <- c(1,3,5,10,50)
cops <- c("Clayton", "Frank", "Gumbel")

dfplot <- data.frame(matrix(ncol = 4, nrow = 0))
names(dfplot) <- c("inputCopula", "fittingMethod", "days", "paramValue")

counterGood <- 0

for (MC_rep in 1:MC_reps) {
  for (copula in cops) {
    for (method in cops) {
      for (n in days) {
        
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
        
        ### Fit copula
        # The input for the copulas are the CDF values of the margins
        
        obs_train_CDF <- c()
        for (i in 1:dim(obs)[2]) {
          obs_train_CDF <- cbind(obs_train_CDF, pnorm(obs[,i], mean = 0, sd = 1))
        }
        
        tryCatch({
        if (method == "Clayton") {
          fitcop <-fitCopula(claytonCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
          params <- fitcop@estimate
          tauValue <- params/(params + 2)
        } else if (method == "Frank") {
          fitcop <-fitCopula(frankCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
          params <- fitcop@estimate
          tauValue <- par2tau("frk", params)
        } else if (method == "Gumbel") {
          fitcop <-fitCopula(gumbelCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
          params <- fitcop@estimate
          tauValue <- par2tau("gum", params)
        } 
        
        
        newdf <- data.frame(inputCopula = copula, fittingMethod = method, days = n, paramValue = params, tauValue = tauValue)
        dfplot <- rbind(dfplot, newdf)
        counterGood <- counterGood + 1
        })
      }
    }
  }
}

sprintf("Success rate is %0.1f%%",100*counterGood/(length(cops)^2*MC_reps*length(days)))

dfplot$daysString <- as.character(dfplot$days)


dfplot$differenceTau <- dfplot$tauValue - tau

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("1" = mypal[1],
               "3" = mypal[2],
               "5" = mypal[3],
              "10" = mypal[4],
              "50" = mypal[5])

  
  

# Plot for tau
quants <- unname(quantile(dfplot$differenceTau, c(0.1, 0.9)))

ylimits <- c(1.5 * quants[1], 1.5 *quants[2])

level_order <- factor(dfplot$daysString, level = as.character(days))


p1 <- ggplot(dfplot, aes(level_order, differenceTau, colour = daysString))
p1 <- p1 + ylim(ylimits[1], ylimits[2])
p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(inputCopula), cols = vars(fittingMethod),
                      labeller = label_bquote(rows = observations: .(inputCopula),
                                              cols = fitting: .(fittingMethod)))



p1 <- p1 + ggtitle(bquote(Evaluation~of~concordance~after~fitting~(tau ==.(tau))))



p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Days") + ylab(bquote(Difference~fitted~and~actual~value~of~tau)) 
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Training days") +
  theme(plot.title = element_text(hjust = 0.5))

plotWidth <- 9
plotHeight <- 6
res <- 400
fileName <- "Arch_ParamFitting_Tau.png"

plot_folder <- paste0("../Data/Plots/")

ggsave(
  paste0(plot_folder, fileName),
  p1,
  width = plotWidth,
  height = plotHeight,
  dpi = res
)


# Plot for theta
dfplot$differenceTheta <- dfplot$paramValue
dfplot$differenceTheta[dfplot$fittingMethod == "Clayton"] <- dfplot$differenceTheta[dfplot$fittingMethod == "Clayton"] - 2 * tau / (1-tau)
dfplot$differenceTheta[dfplot$fittingMethod == "Frank"] <- dfplot$differenceTheta[dfplot$fittingMethod == "Frank"] - tau2par("frk",tau)
dfplot$differenceTheta[dfplot$fittingMethod == "Gumbel"] <- dfplot$differenceTheta[dfplot$fittingMethod == "Gumbel"] - tau2par("gum",tau)

quants <- unname(quantile(dfplot$differenceTheta, c(0.1, 0.9)))

ylimits <- c(1.5 * quants[1], 1.5 *quants[2])

level_order <- factor(dfplot$daysString, level = as.character(days))


p1 <- ggplot(dfplot, aes(level_order, differenceTheta, colour = daysString))
p1 <- p1 + ylim(ylimits[1], ylimits[2])
p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(inputCopula), cols = vars(fittingMethod),
                      labeller = label_bquote(rows = observations: .(inputCopula),
                                              cols = fitting: .(fittingMethod)))



p1 <- p1 + ggtitle(bquote(Evaluation~of~parameter~fitting~(tau==.(tau))))



p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Days") + ylab("Difference fitted parameter and parameter by tau inversion")
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Training days") +
  theme(plot.title = element_text(hjust = 0.5))

plotWidth <- 9
plotHeight <- 6
res <- 400
fileName <- "Arch_ParamFitting_Param.png"

plot_folder <- paste0("../Data/Plots/")

ggsave(
  paste0(plot_folder, fileName),
  p1,
  width = plotWidth,
  height = plotHeight,
  dpi = res
)



