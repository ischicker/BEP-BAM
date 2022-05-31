rm(list=ls())

library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")
groupNR <- 2
fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable


# Save settings
plot_folder <- paste0("../Data/Plots/Group ", groupNR, "/")
dir.create(file.path(plot_folder), showWarnings = FALSE)
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



## PIT diagrams (number of bins = m)

# Get bin number for each row
binNumber <- function(values, obs) {
  sortedValues <- sort(values)
  
  binNum <- 1
  equalBins <- c()
  
  while (binNum <= length(values) & sortedValues[binNum] <= obs) {
    # If multiple hits, select at random
    if (sortedValues[binNum] == obs) {
      equalBins <- c(equalBins, binNum)
    }
    binNum <- binNum + 1
  }
  
  
  if (length(equalBins) > 0) {
    return(sample(equalBins, 1))
  } else {
    return(binNum)
  }
}


# Add list to store bin number
res$bin_list <- list()


mvpp_approaches <- names(res$mvpp_list)
n <- dim(res$mvpp_list$ens)[1]
m <- dim(res$mvpp_list$ens)[2]
d <- dim(res$mvpp_list$ens)[3]

# Add bin number to data
for (model in mvpp_approaches) {
  res$bin_list[[model]] <- array(NA, dim = c(n, d))
  for (day in 1:n) {
    for (i.d in 1:d) {
      if (model %in% c("clayton", "frank", "gumbel")){
        dat <- res$mvpp_list[[model]]$mvppout
        if (!any(is.na(dat[day,,i.d]))) {
          res$bin_list[[model]][day, i.d] <- binNumber(dat[day,,i.d], res$obs[day, i.d])
        }
      } else {
        res$bin_list[[model]][day, i.d] <- binNumber(res$mvpp_list[[model]][day,,i.d], res$obs[day, i.d])
      }
    }
  }
}


# Create the PIT plots
createPIT <- function(model, i.d) {
  dat <- res$bin_list[[model]][,i.d]
  nonNANumber <- sum(!is.na(dat))
  height <- nonNANumber / (m + 1)
  highestValue <- max(sapply(1:(m+1), FUN = function(x) length(dat[dat == x])))
  p <- ggplot() + aes(dat) + geom_histogram(breaks = seq(0, m + 1, 1)) +
    geom_hline(yintercept = height, col="steelblue", linetype = "dashed") + 
    scale_x_continuous(breaks = seq(0, m + 1, 1), labels=round(seq(0, m + 1, 1)/ (m + 1), 2)) + 
    scale_y_continuous(breaks = seq(0, highestValue, height), labels = seq(0, highestValue / height, 1)) +
    labs(y = "Frequency ratio", x = "PIT value") +
    ggtitle(paste0(model," for d = ", i.d)) + theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

cols <- 4
rows <- 3



# Create and save the PIT plots
for (i.d in 1:d) {
  plot_vec <- c()
  plotWidth <-  8
  plotHeight <- 8
  for (model in mvpp_approaches) {
    plot <- createPIT(model, i.d)
    plot_vec <- c(plot_vec, list(plot))
    # Save the plots
    savePlots(paste0("PIT_group_", groupNR, "_", model,"_d_", i.d, ".png"), plot)
  }
  plotWidth <-  5 * cols
  plotHeight <- 5 * rows
  savePlots(paste0("PIT_group_", groupNR, "_grid_d_",i.d,".png"),ggarrange(plotlist = plot_vec,nrow = rows,ncol = cols))
}





