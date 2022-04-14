require(MASS)
require(copula)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")


MC_reps <- 100
theta0 <- 5
d <- 3
days <- c(1,3,5,10,50)
cops <- c("Clayton", "Frank", "Gumbel")

dfplot <- data.frame(matrix(ncol = 4, nrow = 0))
names(dfplot) <- c("inputCopula", "fittingMethod", "days", "paramValue")

for (MC_rep in 1:MC_reps) {
  for (copula in cops) {
    for (method in cops) {
      for (n in days) {
        
        # Generate observations
        if (copula == "Frank"){
          # -infty < theta0 < infty and theta0 != 0
          cop <- frankCopula(param = theta0, dim = d)
        } else if (copula == "Gumbel") {
          # 1 <= theta0 < infty
          cop <- gumbelCopula(param = theta0, dim = d)
        } else if (copula == "Clayton") {
          # -1 <= theta0 < infty and theta0 != 0
          cop <- claytonCopula(param = theta0, dim = d)
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
        
        try({
        if (method == "Clayton") {
          fitcop <-fitCopula(claytonCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        } else if (method == "Frank") {
          fitcop <-fitCopula(frankCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        } else if (method == "Gumbel") {
          fitcop <-fitCopula(gumbelCopula(dim = d), data = obs_train_CDF, method="mpl", start = 1,optim.control = list(maxit=1000), upper = 100)
        } 
        params <- fitcop@estimate
        newdf <- data.frame(inputCopula = copula, fittingMethod = method, days = n, paramValue = params)
        dfplot <- rbind(dfplot, newdf)
        })
      }
    }
  }
}

dfplot$daysString <- as.character(dfplot$days)

dfplot$difference <- dfplot$paramValue - theta0

mypal <- colorspace::rainbow_hcl(5)
mypal_use <- c("1" = mypal[1],
               "3" = mypal[2],
               "5" = mypal[3],
              "10" = mypal[4],
              "50" = mypal[5])

  
  
  
quants <- unname(quantile(dfplot$difference, c(0.1, 0.9)))

ylimits <- c(1.5 * quants[1], 1.5 *quants[2])

level_order <- factor(dfplot$daysString, level = as.character(days))


p1 <- ggplot(dfplot, aes(level_order, difference, colour = daysString))
p1 <- p1 + ylim(ylimits[1], ylimits[2])
p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
p1 <- p1 + facet_grid(rows = vars(inputCopula), cols = vars(fittingMethod),
                      labeller = label_bquote(rows = observations: .(inputCopula),
                                              cols = fitting: .(fittingMethod)))



p1 <- p1 + ggtitle(bquote(Evaluation~of~parameter~fitting~(theta[0]==.(theta0))))



p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Days") + ylab("Difference fitted and actual parameter") 
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Training days") +
  theme(plot.title = element_text(hjust = 0.5))

plotWidth <- 9
plotHeight <- 6
res <- 400
fileName <- "Arch_ParamFitting.png"

plot_folder <- paste0("../Data/Plots/")

ggsave(
  paste0(plot_folder, fileName),
  p1,
  width = plotWidth,
  height = plotHeight,
  dpi = res
)






