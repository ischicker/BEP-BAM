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
# "source" directory
setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/simulation code")
source("getData.R")
getData()


library(ggplot2)
library(ggpubr)
library(purrr)
library(qqplotr)

plotTheme <- theme(axis.text=element_text(size=20),
                   axis.title=element_text(size=25,face="bold"),
                   plot.title = element_text(size=35,hjust = 0.5))

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
res <- 250

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
    map(c("obs", ensembleMembers), function(member) {
      ggplot(data, aes_string(x = "time", y = member)) + geom_point() +
        labs(x = "Time", y = "Temperature (K)") +
        ggtitle(member) + plotTheme
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
  observationsPlot <- list(ggplot(data, aes(x=obs)) + geom_histogram())
  
  return(c(observationsPlot, ensemblePlots))
}

# Histograms of all values
histPlots1 <- createHistPlots(data1)
histPlots2 <- createHistPlots(data2)
histPlots3 <- createHistPlots(data3)

# Data from a single month
data1month <- data1[1:30,]
data2month <- data2[1:30,]
data3month <- data3[1:30,]

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
    annotate("text", size = 12, x=max(data[[member]]),y=min(data[[member]]),hjust=1,label=paste0("Shapiro p-value: ", signif(shapiro.test(data[[member]])$p.value, digits = 3))) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggtitle(member) + plotTheme
})
  nonNAdata <- data$obs[!is.na(data$obs)]
  observationsPlot <- list(ggplot(data = data, mapping = aes(sample = obs)) +
    stat_qq_band() +
    stat_qq_line() +
    stat_qq_point() +
    annotate("text",size = 12,x=max(nonNAdata),y=min(nonNAdata),hjust=1,label=paste0("Shapiro p-value: ", signif(shapiro.test(nonNAdata)$p.value, digits = 3))) +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggtitle("obs") + plotTheme )
  
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

savePlots("SeasonPair_group_2.png",ggarrange(plotlist = seasonPlot2[1:2],nrow = 1,ncol = 2))

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

savePlots("QQPair_group_2.png",ggarrange(plotlist = rbind(qqPlots2[1],qqPlots2month[1]),nrow = 1,ncol = 2))

## Investigate correlation between observations and forecasts
createCorrelationPlots <- function(data, compareTo) {
  
  return(
    map(c("obs", ensembleMembers), function(member) {
      ggplot(data, aes_string(x = compareTo, y = member)) + geom_point() +
        geom_abline(intercept = 0, slope = 1, color="steelblue", 
                    linetype="dashed", size = 1.5)  + plotTheme
        
    })
  )
}

# Get the plots
correlationPlots1 <- createCorrelationPlots(data1, "obs")
correlationPlots2 <- createCorrelationPlots(data2, "obs")
correlationPlots3 <- createCorrelationPlots(data3, "obs")

# Correlations between ensemble members
correlationPlotsEns1 <- createCorrelationPlots(data1, "laef1")
correlationPlotsEns2 <- createCorrelationPlots(data2, "laef1")
correlationPlotsEns3 <- createCorrelationPlots(data3, "laef1")

# Save the plots
savePlots("CorrelationGrid_group_1.png",ggarrange(plotlist = correlationPlots1,nrow = rows,ncol = cols))
savePlots("CorrelationGrid_group_2.png",ggarrange(plotlist = correlationPlots2,nrow = rows,ncol = cols))
savePlots("CorrelationGrid_group_3.png",ggarrange(plotlist = correlationPlots3,nrow = rows,ncol = cols))

savePlots("CorrelationPair_group_2.png",ggarrange(plotlist = correlationPlots2[1:2],nrow = 1,ncol = 2))
savePlots("CorrelationPane_group_2.png",ggarrange(plotlist = correlationPlots2[2],nrow = 1,ncol = 1))

savePlots("CorrelationGridEnsemble_group_1.png",ggarrange(plotlist = correlationPlotsEns1,nrow = rows,ncol = cols))
savePlots("CorrelationGridEnsemble_group_2.png",ggarrange(plotlist = correlationPlotsEns2,nrow = rows,ncol = cols))
savePlots("CorrelationGridEnsemble_group_3.png",ggarrange(plotlist = correlationPlotsEns3,nrow = rows,ncol = cols))

savePlots("CorrelationPairEnsemble_group_2.png",ggarrange(plotlist = correlationPlotsEns2[2:3],nrow = 1,ncol = 2))
savePlots("CorrelationPaneEnsemble_group_2.png",ggarrange(plotlist = correlationPlotsEns2[3],nrow = 1,ncol = 1))

savePlots("CorrelationPanes.png",ggarrange(plotlist = c(correlationPlots2[2],correlationPlotsEns2[3]),nrow = 1,ncol = 2))

plotWidth <- 10
plotHeight <- 10
savePlots("CorrelationPane_group_2.png",ggarrange(plotlist = correlationPlots2[2],nrow = 1,ncol = 1))
savePlots("CorrelationPaneEnsemble_group_2.png",ggarrange(plotlist = correlationPlotsEns2[3],nrow = 1,ncol = 1))
