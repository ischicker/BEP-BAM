rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")


source("../Settings.R")

getModelSettings(modelSetting = 1)

setting <- 1

load(paste0("../Data/TestStatistic/TestStatistic_Archimedean","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".Rdata"))

df <- dfmc; rm(dfmc)

df1 <- subset(df)


df1$value <- (-1)*df1$value
df0 <- df1

input_scores <- unique(df1$score)
plot_folder <- paste0("../Data/Plots/Arch","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,"/")
dir.create(file.path(plot_folder), showWarnings = FALSE)


# drop model == "ens"
df1_save <- df1 
df1 <- subset(df1_save, model != "ens")

# also drop EMOS.Q
df2 <- subset(df1, model != "emos.q")

mypal <- colorspace::rainbow_hcl(9)
mypal_use <- c("decc.q" = mypal[1],
               "ecc.q" = mypal[2],
               "ecc.s" = mypal[3],
               "gca" = mypal[4],
               "ssh" = mypal[5],
               "Clayton" = mypal[6],
               "Frank" = mypal[7],
               "Gumbel" = mypal[8],
               "GOF" = mypal[9])

df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh", "Clayton", "Frank", "Gumbel", "GOF"))
model_vec <- c("dECC", "ECC-Q", "ECC-S", "GCA", "SSh", "Clayton", "Frank", "Gumbel", "GOF")

ylimitFunc <- function(val1, val2) {
  return(1.5 * max(abs(val1), abs(val2)))
}


# Plot the scores for copula
plotScores <- function(dfplot, cop, this_score){
  

  dfplotCop <- subset(dfplot, copula == cop)
  
  alpha <- 0.25
  
  quants <- unname(quantile(dfplot$value, c(0.01, 0.99)))
  
  ylimits <- c(1.5 * min(quants[1], qnorm(alpha)), 1.5 * max(quants[2], qnorm(1 - alpha)))
  

  
  p1 <- ggplot(dfplotCop, aes(model, value, colour = model))
  p1 <- p1 + ylim(ylimits[1], ylimits[2])
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(alpha), ymax=qnorm(1-alpha)), fill = "gray75", color="gray75", alpha=alpha)
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  p1 <- p1 + facet_grid(rows = vars(theta), cols = vars(theta0),
                        labeller = label_bquote(rows = theta==.(theta),
                                                cols = theta[0]==.(theta0)))
  
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "es"){ 
    title <- paste("Energy Score for",cop,"copula")
    p1 <- p1 + ggtitle(title)
  } 
  if(scval == "vs1"){
    title <- paste("Variogram Score for",cop,"copula")
    p1 <- p1 + ggtitle(title)
  }
  if(scval == "vs0"){
    title <- paste("Variogram Score (p=0.5) for",cop,"copula")
    p1 <- p1 + ggtitle(title)
  }
  
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  return(p1)
}

plotScoresTau <- function(dfplot, cop, this_score, isAverage){
  
  
  dfplotCop <- subset(dfplot, copula == cop)
  
  alpha <- 0.25
  
  quants <- unname(quantile(dfplot$value, c(0.01, 0.99)))
  
  ylimits <- c(1.5 * min(quants[1], qnorm(alpha)), 1.5 * max(quants[2], qnorm(1 - alpha)))
  
  
  
  p1 <- ggplot(dfplotCop, aes(model, value, colour = model))
  p1 <- p1 + ylim(ylimits[1], ylimits[2])
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(alpha), ymax=qnorm(1-alpha)), fill = "gray75", color="gray75", alpha=alpha)
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  
  if (!isAverage){
    p1 <- p1 + facet_grid(rows = vars(tau),
                          labeller = label_bquote(rows = tau==.(tau)))
  }
  
  scval <- strsplit(this_score, split = "_")[[1]][1]
  print(scval)
  if(scval == "es"){ 
    title <- paste("Energy Score for",cop,"copula")
    p1 <- p1 + ggtitle(title)
  } 
  if(scval == "vs1"){
    title <- paste("Variogram Score (p=1) for",cop,"copula")
    p1 <- p1 + ggtitle(title)
  }
  if(scval == "vs0"){
    title <- paste("Variogram Score (p=0.5) for",cop,"copula")
    # print(title)
    p1 <- p1 + ggtitle(title)
  }
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  return(p1)
}

plotScoresLerch <- function(dfplot, this_score){
  
  
  alpha <- 0.25
  
  quants <- unname(quantile(dfplot$value, c(0.01, 0.99)))
  
  ylimits <- c(1.5 * min(quants[1], qnorm(alpha)), 1.5 * max(quants[2], qnorm(1 - alpha)))
  
  
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model))
  p1 <- p1 + ylim(ylimits[1], ylimits[2])
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(alpha), ymax=qnorm(1-alpha)), fill = "gray75", color="gray75", alpha=alpha)
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  
  
  p1 <- p1 + facet_grid(rows = vars(rho), cols = vars(rho0),
                        labeller = label_bquote(rows = rho==.(rho),
                                                cols = rho[0] == .(rho0)))

  
  scval <- strsplit(this_score, split = "_")[[1]][1]
  if(scval == "es"){ 
    title <- bquote("Energy Score,"~sigma==sqrt(5))
    p1 <- p1 + ggtitle(title)
  } 
  if(scval == "vs1"){
    title <- bquote("Variogram Score (p=1),"~sigma==sqrt(5))
    p1 <- p1 + ggtitle(title)
  }
  if(scval == "vs0"){
    title <- bquote("Variogram Score (p=0.5),"~sigma==sqrt(5))
    p1 <- p1 + ggtitle(title)
  }
  
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
  p1 <- p1 + scale_x_discrete(label = model_vec)
  return(p1)
}

saveFigure <- function(fileName, fig) {
  res <- 250
  
  ggsave(
    paste0(plot_folder, fileName),
    fig,
    width = plotWidth,
    height = plotHeight,
    dpi = res,
    limitsize = FALSE
  )
}

## ES
this_score <- "es_list"
dfplot <- subset(df2, score == this_score)

if (observationsModel == 1) {
  
  plotWidth <- 9
  plotHeight <- 6

  dfplot$equalThetas <- (dfplot$theta0 == dfplot$theta)
  # Get the plots
  p1saveFrank <- plotScores(dfplot, "Frank", this_score)
  p1saveClayton <- plotScores(dfplot, "Clayton", this_score)
  p1saveGumbel <- plotScores(dfplot, "Gumbel", this_score)
  
  # Save the plots individually
  saveFigure(paste0("Arch_ES_Clayton","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveClayton)
  saveFigure(paste0("Arch_ES_Frank","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveFrank)
  saveFigure(paste0("Arch_ES_Gumbel","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveGumbel)
  print("Finished individual plots ES")
} else if (observationsModel == 2) {
  
  # plotWidth <- 9
  # plotHeight <- 6 * max(dfplot$repetition)
  # 
  # p1saveFrank <- plotScoresTau(dfplot, "Frank", this_score, FALSE)
  # p1saveClayton <- plotScoresTau(dfplot, "Clayton", this_score, FALSE)
  # p1saveGumbel <- plotScoresTau(dfplot, "Gumbel", this_score,FALSE)
  # 
  # # Save the plots individually
  # saveFigure(paste0("Arch_ES_Separated_Clayton","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveClayton)
  # saveFigure(paste0("Arch_ES_Separated_Frank","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveFrank)
  # saveFigure(paste0("Arch_ES_Separated_Gumbel","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveGumbel)
  # 
  # print("Finished ES Separated")
  
  plotWidth <- 9
  plotHeight <- 6
  
  p1saveFrank <- plotScoresTau(dfplot, "Frank", this_score, TRUE)
  p1saveClayton <- plotScoresTau(dfplot, "Clayton", this_score, TRUE)
  p1saveGumbel <- plotScoresTau(dfplot, "Gumbel", this_score,TRUE)
  
  # Save the plots individually
  saveFigure(paste0("Arch_ES_Together_Clayton","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveClayton)
  saveFigure(paste0("Arch_ES_Together_Frank","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveFrank)
  saveFigure(paste0("Arch_ES_Together_Gumbel","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveGumbel)
  
  print("Finished ES Together")
  
} else if (observationsModel == 3) {
  plotWidth <- 12
  plotHeight <- 6 
  
  p1 <- plotScoresLerch(dfplot, this_score)
  saveFigure(paste0("Lerch_ES","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1)
}



# if (observationsModel == 1) {
#   ### Join the plots
#   
#   library(gridExtra)
#   plotWidth <- 9
#   plotHeight <- 25
#   saveFigure(paste0("Arch_ES_Combined","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"), grid.arrange(p1saveClayton,p1saveFrank,p1saveGumbel, ncol = 1))
#   print("Finished combined plots ES")
# }
## VS

this_score <- "vs0_list"

dfplot <- subset(df2, score == this_score)


if (observationsModel == 1) {
  
  plotWidth <- 9
  plotHeight <- 6
  
  dfplot$equalThetas <- (dfplot$theta0 == dfplot$theta)
  # Get the plots
  p1saveFrank <- plotScores(dfplot, "Frank", this_score)
  p1saveClayton <- plotScores(dfplot, "Clayton", this_score)
  p1saveGumbel <- plotScores(dfplot, "Gumbel", this_score)
  
  # Save the plots individually
  saveFigure(paste0("Arch_VS_Clayton","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveClayton)
  saveFigure(paste0("Arch_VS_Frank","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveFrank)
  saveFigure(paste0("Arch_VS_Gumbel","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveGumbel)
  print("Finished individual plots VS")
} else if (observationsModel == 2) {
  
  # plotWidth <- 9
  # plotHeight <- 6 * max(dfplot$repetition)
  # 
  # p1saveFrank <- plotScoresTau(dfplot, "Frank", this_score, FALSE)
  # p1saveClayton <- plotScoresTau(dfplot, "Clayton", this_score, FALSE)
  # p1saveGumbel <- plotScoresTau(dfplot, "Gumbel", this_score,FALSE)
  # 
  # # Save the plots individually
  # saveFigure(paste0("Arch_VS_Separated_Clayton","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveClayton)
  # saveFigure(paste0("Arch_VS_Separated_Frank","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveFrank)
  # saveFigure(paste0("Arch_VS_Separated_Gumbel","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveGumbel)
  # 
  # print("Finished VS Separated")
  
  plotWidth <- 9
  plotHeight <- 6
  
  p1saveFrank <- plotScoresTau(dfplot, "Frank", this_score, TRUE)
  p1saveClayton <- plotScoresTau(dfplot, "Clayton", this_score, TRUE)
  p1saveGumbel <- plotScoresTau(dfplot, "Gumbel", this_score,TRUE)
  
  if (this_score == "vs0_list") {
    text <- "_weight_0_p_05_"
  }
  
  # Save the plots individually
  saveFigure(paste0("Arch_VS", text, "Together_Clayton","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveClayton)
  saveFigure(paste0("Arch_VS", text, "Together_Frank","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveFrank)
  saveFigure(paste0("Arch_VS", text, "Together_Gumbel","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1saveGumbel)
  
  print("Finished VS Together")
  
} else if (observationsModel == 3) {
  plotWidth <- 12
  plotHeight <- 6 
  
  p1 <- plotScoresLerch(dfplot, this_score)
  saveFigure(paste0("Lerch_VS","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"),p1)
}

# if (observationsModel == 1) {
#   ### Join the plots
#   
#   plotWidth <- 9
#   plotHeight <- 25 
#   saveFigure(paste0("Arch_VS_Combined","_setting_",setting, "_obsmodel_",observationsModel,"_fcmodel_",forecastModel,".png"), grid.arrange(p1saveClayton,p1saveFrank,p1saveGumbel, ncol = 1))
#   print("Finished combined plots VS")
# }