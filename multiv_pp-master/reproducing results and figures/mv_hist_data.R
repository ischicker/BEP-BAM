rm(list=ls())
library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")
groupNR <- 2
fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable

dat <- res$mvpp_list

days <- dim(dat$ens)[1]
m <- dim(dat$ens)[2]

# Save settings
# Create and save the Multivariate PIT plots
plot_folder <- paste0("../Data/Plots/Group ", groupNR, "/Multivariate/")
dir.create(file.path(plot_folder), showWarnings = FALSE)
plot_vec <- c()
plotWidth <-  8
plotHeight <- 8
resolution <- 400

# Function to save the plots
savePlots <- function(fileName, plot) {
  ggsave(
    paste0(plot_folder, fileName),
    plot,
    width = plotWidth,
    height = plotHeight,
    dpi = resolution
  )
}

reformat <- function(x, modelName)
{
  B <- list()
  
  if (modelName %in% c("clayton", "frank", "gumbel")) {
    x <- x$mvppout
  }
  for(i in 1:days){
    B[[i]] <- matrix(x[i,,], nrow = dim(x)[2], ncol = dim(x)[3])

  }
  
  return(B)
}

## Multivariate ranks 
mv.rank <- function(x)
{
  d <- dim(x)
  x.prerank <- numeric(d[2])
  for(i in 1:d[2]) {
    x.prerank[i] <- sum(apply(x<=x[,i],2,all))
  }
  x.rank <- rank(x.prerank,ties="random")
  return(x.rank)
}


## Multivariate rank histograms
mvr.histogram <- function(B,modelName)
{
  x <- c() 
  for(i in 1:days){
    x <- c(x, mv.rank(t(B[[i]]))[1])
  }
  
  # .GlobalEnv$x <- x
  
  # hist(x,breaks=seq(0,20,by=1),main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray40",border="white",ylim=hist_ylim)
  # ggplot(dfplot, aes(model, value, colour = model))
  intercept <- days/m
  highestValue <- max(sapply(1:m, FUN = function(s) sum(x == s)))

  p <- ggplot() + aes(x) + geom_histogram(breaks = seq(0, m, 1)) +
    geom_hline(yintercept = intercept, col="steelblue", linetype = "dashed") + 
    scale_x_continuous(breaks = seq(0, m , 1), labels=round(seq(0, m , 1)/ (m), 2)) + 
    scale_y_continuous(breaks = seq(0, highestValue, intercept), labels = seq(0, highestValue / intercept, 1)) +
    labs(y = "Frequency ratio", x = "Normalized rank value")
  

  p <- p + ggtitle(paste("Multivariate rank histogram for", modelName)) + theme(plot.title = element_text(hjust = 0.5))
    
  return(p)

}

plotWidth <-  8
plotHeight <- 8
plot_vec <- c()
for (model in names(dat)) {
  if (model != "emos.q") {
    y <- reformat(dat[[model]], model)
    p <- mvr.histogram(y,modelName = model)
    plot_vec <- c(plot_vec, list(p))
    savePlots(paste0("PIT_group_", groupNR, "_", model,"_multivariate.png"), p)
  }
}

cols <- 3
rows <- 3


plotWidth <-  5 * cols
plotHeight <- 5 * rows
savePlots(paste0("PIT_group_", groupNR, "_grid_multivariate.png"),ggarrange(plotlist = plot_vec,nrow = rows,ncol = cols))




