### Input consists of
# 72 different lead times
# 216 stations (input$leadtimen consists of the stations)
# 17 ensemble members (laef1, ..., leaf17)
# input$leadtimen$lt contains the lead time n
# 990 different time stamps (unique values of "initial")
# td contains the discretized time values
# initial is the time of measurement
# ytime = initial + 1
# inca (Integrated Nowcasting through Comprehensive Analyses) contains the observations
# cosmo (Consortium for Small-scale Modeling) another forecast

rm(list=ls())

library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)


# "source" directory 
setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/simulation code")
source("getData.R")
getData()




## Create a grid of histograms of ensemble members

# Number of ensemble members
m <- sum(grepl("laef", names(data1)))

# Plots are put in a grid
cols <- 6
rows <- 3

ensembleMembers <- sapply(1:m, FUN = function(x) paste0("laef", x))


# Save settings
saveFolder <- "Combined_LAEF_Data"
plot_folder <- paste0("../Data/Plots/", saveFolder, "/")
dir.create(file.path(plot_folder), showWarnings = FALSE)
dir.create(file.path(paste0(plot_folder, "Month/")), showWarnings = FALSE)
plotWidth <- cols * 4
plotHeight <- rows * 4
res <- 400

# Function to save the plots
savePlots <- function(fileName, plot) {
  ggsave(
    paste0(plot_folder, fileName),
    plot,
    width = plotWidth,
    height = plotHeight,
    dpi = res
  )
}

## Investigate temporal dependence of data for one particular station
seasonData <- subset(data1, stat %in% statNR[group1[1]])

# Create plots for T over time
createSeasonPlots <- function(data) {
  
  return(
    map(c("inca", ensembleMembers), function(member) {
      ggplot(data, aes_string(x = "time", y = member)) + geom_point() +
        labs(x = "Time", y = "Temperature (K)") +
        ggtitle(member) + theme(plot.title = element_text(hjust = 0.5))
    })
  )
}

# Get the plots
seasonPlot1 <- createSeasonPlots(data1)
seasonPlot2 <- createSeasonPlots(data2)
seasonPlot3 <- createSeasonPlots(data3)


# Create histograms
createHistPlots <- function(data){
  ensemblePlots <- map(ensembleMembers, function(member) {ggplot(data, aes_string(x=member)) + geom_histogram()})
  observationsPlot <- list(ggplot(data, aes(x=inca)) + geom_histogram())
  
  return(c(observationsPlot, ensemblePlots))
}

# Histograms of all values
histPlots1 <- createHistPlots(data1)
histPlots2 <- createHistPlots(data2)
histPlots3 <- createHistPlots(data3)

# Histograms of one month
histPlots1month <- createHistPlots(data1month)
histPlots2month <- createHistPlots(data2month)
histPlots3month <- createHistPlots(data3month)

# QQ plots on normality assumption with 95% confidence interval
createQqPlots <- function(data){
  
  ensemblePlots <- map(ensembleMembers, function(member) {
  ggplot(data = data, mapping = aes_string(sample = member)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    annotate("text",x=max(data[[member]]),y=min(data[[member]]),hjust=1,label=paste0("Shapiro p-value: ", signif(shapiro.test(data[[member]])$p.value, digits = 3))) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggtitle(member) + theme(plot.title = element_text(hjust = 0.5))
})
  nonNAdata <- data$inca[!is.na(data$inca)]
  observationsPlot <- list(ggplot(data = data, mapping = aes(sample = inca)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    annotate("text",x=max(nonNAdata),y=min(nonNAdata),hjust=1,label=paste0("Shapiro p-value: ", signif(shapiro.test(nonNAdata)$p.value, digits = 3))) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggtitle("inca") + theme(plot.title = element_text(hjust = 0.5)) )
  
  return(c(observationsPlot, ensemblePlots))
}

# QQ plots of all data
qqPlots1 <- createQqPlots(data1)
qqPlots2 <- createQqPlots(data2)
qqPlots3 <- createQqPlots(data3)

# QQ plots of one month
qqPlots1month <- createQqPlots(data1month)
qqPlots2month <- createQqPlots(data2month)
qqPlots3month <- createQqPlots(data3month)


# Save the plots
savePlots("SeasonGrid_group_1.png",ggarrange(plotlist = seasonPlot1,nrow = rows,ncol = cols))
savePlots("SeasonGrid_group_2.png",ggarrange(plotlist = seasonPlot2,nrow = rows,ncol = cols))
savePlots("SeasonGrid_group_3.png",ggarrange(plotlist = seasonPlot3,nrow = rows,ncol = cols))

# Save histograms
savePlots("HistogramGrid_group_1.png",ggarrange(plotlist = histPlots1,nrow = rows,ncol = cols))
savePlots("HistogramGrid_group_2.png",ggarrange(plotlist = histPlots2,nrow = rows,ncol = cols))
savePlots("HistogramGrid_group_3.png",ggarrange(plotlist = histPlots3,nrow = rows,ncol = cols))

savePlots("Month/HistogramGrid_group_1.png",ggarrange(plotlist = histPlots1month,nrow = rows,ncol = cols))
savePlots("Month/HistogramGrid_group_2.png",ggarrange(plotlist = histPlots2month,nrow = rows,ncol = cols))
savePlots("Month/HistogramGrid_group_3.png",ggarrange(plotlist = histPlots3month,nrow = rows,ncol = cols))


# Save QQplots
savePlots("QQGrid_group_1.png",ggarrange(plotlist = qqPlots1,nrow = rows,ncol = cols))
savePlots("QQGrid_group_2.png",ggarrange(plotlist = qqPlots2,nrow = rows,ncol = cols))
savePlots("QQGrid_group_3.png",ggarrange(plotlist = qqPlots3,nrow = rows,ncol = cols))

savePlots("Month/QQGrid_group_1.png",ggarrange(plotlist = qqPlots1month,nrow = rows,ncol = cols))
savePlots("Month/QQGrid_group_2.png",ggarrange(plotlist = qqPlots2month,nrow = rows,ncol = cols))
savePlots("Month/QQGrid_group_3.png",ggarrange(plotlist = qqPlots3month,nrow = rows,ncol = cols))

## Investigate correlation between observations and forecasts
createCorrelationPlots <- function(data, compareTo) {
  
  return(
    map(c("inca", ensembleMembers), function(member) {
      ggplot(data, aes_string(x = compareTo, y = member)) + geom_point() +
        geom_abline(intercept = 0, slope = 1, color="steelblue", 
                    linetype="dashed", size = 1.5)
    })
  )
}

# Get the plots
correlationPlots1 <- createCorrelationPlots(data1, "inca")
correlationPlots2 <- createCorrelationPlots(data2, "inca")
correlationPlots3 <- createCorrelationPlots(data3, "inca")

# Correlations between ensemble members
correlationPlotsEns1 <- createCorrelationPlots(data1, "laef1")
correlationPlotsEns2 <- createCorrelationPlots(data2, "laef1")
correlationPlotsEns3 <- createCorrelationPlots(data3, "laef1")

# Save the plots
savePlots("CorrelationGrid_group_1.png",ggarrange(plotlist = correlationPlots1,nrow = rows,ncol = cols))
savePlots("CorrelationGrid_group_2.png",ggarrange(plotlist = correlationPlots2,nrow = rows,ncol = cols))
savePlots("CorrelationGrid_group_3.png",ggarrange(plotlist = correlationPlots3,nrow = rows,ncol = cols))

savePlots("CorrelationGridEnsemble_group_1.png",ggarrange(plotlist = correlationPlotsEns1,nrow = rows,ncol = cols))
savePlots("CorrelationGridEnsemble_group_2.png",ggarrange(plotlist = correlationPlotsEns2,nrow = rows,ncol = cols))
savePlots("CorrelationGridEnsemble_group_3.png",ggarrange(plotlist = correlationPlotsEns3,nrow = rows,ncol = cols))


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

# Get list of bin numbers divided by number
binVector <- function(dat) {
  return(
    apply(dat[c(ensembleMembers, "obs")], MARGIN = 1, FUN = function(x) {
    if (!any(is.na(x))) {
      ensembleVector <- x[1:m]
      return(binNumber(ensembleVector, x[m+1]))
    } else {
      return(NA)
    }
  }))
}

# Add bin number to data
data1$bin <- binVector(data1)
data2$bin <- binVector(data2)
data3$bin <- binVector(data3)


# Create the PIT plots
createPIT <- function(dat) {
  nonNANumber <- sum(!is.na(dat$bin))
  height <- nonNANumber / (m + 1)
  highestValue <- max(sapply(1:(m+1), FUN = function(x) length(dat$bin[dat$bin == x])))
  p <- ggplot(dat, aes(x=bin)) + geom_histogram(breaks = seq(0, m + 1, 1)) +
    geom_hline(yintercept = height, col="steelblue", linetype = "dashed") + 
    scale_x_continuous(breaks = seq(0, m + 1, 1), labels=round(seq(0, m + 1, 1)/ (m + 1), 2)) + 
    scale_y_continuous(breaks = seq(0, highestValue, height), labels = seq(0, highestValue / height, 1)) +
    labs(y = "Frequency ratio", x = "PIT value")
  
  return(p)
}

# Save the plots
savePlots("PIT_group_1.png",createPIT(data1))
savePlots("PIT_group_3.png",createPIT(data3))

stations1 <- unique(data1$stat)
for (s in 1:length(stations1)) {
  print(s)
  savePlots(paste0("PIT_group_1_d_",s,".png"),createPIT(subset(data1, stat == stations1[s])))
}

stations2 <- unique(data2$stat)
for (s in 1:length(stations2)) {
  print(s)
  savePlots(paste0("PIT_group_2_d_",s,".png"),createPIT(subset(data2, stat == stations2[s])))
}

stations3 <- unique(data3$stat)
for (s in 1:length(stations3)) {
  print(s)
  savePlots(paste0("PIT_group_3_d_",s,".png"),createPIT(subset(data3, stat == stations3[s])))
}
