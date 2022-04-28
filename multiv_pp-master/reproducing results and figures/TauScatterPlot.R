rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")

source("../simulation code/sourceArchimedean/CopulaParameter.R")
source( "../Settings.R")

setting <- 1

getModelSettings(modelSetting = 2)

fName <- paste0("Archimedean","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,"_ID_")


flist <- list.files("../Data/Rdata/")
existing <- as.numeric(sapply(flist, FUN = function(x) as.numeric(strsplit(strsplit(x, fName)[[1]][2], ".Rdata"))))


load(paste0("../Data/Rdata/",fName, "1.Rdata"))

df <- res; rm(res)

df <- list(frankFit = df$param_list$frank, gumbelFit = df$param_list$gumbel, claytonFit = df$param_list$clayton, tau = df$tau)

# df$frankTau <- sapply(df$frankFit, FUN = paramToTau, copula = "Frank")
# df$gumbelTau <- sapply(df$gumbelFit, FUN = paramToTau, copula = "Gumbel")
# df$claytonTau <- sapply(df$claytonFit, FUN = paramToTau, copula = "Clayton")
# 
# dfplot <- data.frame(claytonTau = df$claytonTau, frankTau = df$frankTau, gumbelTau = df$gumbelTau, tau = rep(df$tau, length(df$frankTau)))

dfplot <- data.frame()

for(ID in existing[!is.na(existing)]){
  load(paste0("../Data/Rdata/", fName, ID,".Rdata"))
  res$frankTau <- sapply(res$param_list$frank, FUN = paramToTau, copula = "Frank")
  res$gumbelTau <- sapply(res$param_list$gumbel, FUN = paramToTau, copula = "Gumbel")
  res$claytonTau <- sapply(res$param_list$clayton, FUN = paramToTau, copula = "Clayton") 
  
  new_df <- data.frame(claytonTau = res$claytonTau, frankTau = res$frankTau, gumbelTau = res$gumbelTau, 
                       tau = rep(res$tau, length(res$frankTau)), input_cop = input_par$copula[ID])
  
  dfplot <- rbind(dfplot, new_df)
}

tauValues <- unique(dfplot$tau)


plotTau <- function(fitCopula) {
  p1 <- ggplot(dfplot, aes_string(x="tau", y=fitCopula))
  
  for (i in 1:length(tauValues)) {
    p1 <- p1 +geom_boxplot(data=subset(dfplot, tau == tauValues[i]), width=0.5/length(tauValues))
  }
  p1 <- p1 + geom_abline(intercept = 0, slope = 1, color="grey", 
                    linetype="dashed", size=0.5)
  p1 <- p1 + facet_grid(rows = vars(input_cop),
                        labeller = label_bquote(rows = observations: .(as.character(input_cop))))
  p1 <- p1 + scale_x_continuous(breaks=seq(0,1,0.1)) + scale_y_continuous(breaks=seq(0,1,0.2))
  
  
  p1 <- p1 + theme_bw()
  p1 <- p1 + xlab(bquote(tau[input])) + ylab(bquote(tau[output]))
  
  if (fitCopula == "claytonTau") {
    copula <- "Clayton"
  } else if (fitCopula == "frankTau") {
    copula <- "Frank"
  } else if (fitCopula == "gumbelTau") {
    copula <- "Gumbel"
  } else {
    stop("No copula")
  }
  
  p1 <- p1 + ggtitle(bquote("Input versus fitted" ~ tau ~ "for" ~ .(copula) ~ " fitting")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p1)
}


plot_folder <- paste0("../Data/Plots/Arch_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,"/")
dir.create(file.path(plot_folder), showWarnings = FALSE)

plotWidth <- 9
plotHeight <- 6
res <- 400

for (fitCopula in c("claytonTau", "frankTau", "gumbelTau")){
  fileName <- paste0("Arch_",fitCopula,"_setting_", setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png")
  
  p1 <- plotTau(fitCopula)
  
  ggsave(
    paste0(plot_folder, fileName),
    p1,
    width = plotWidth,
    height = plotHeight,
    dpi = res
  )
}


