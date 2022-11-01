# rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/University/Bachelor/Year 3/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")


groupNR <- 3
fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable

load(paste0("../Data/TestStatistic/TestStatistic_data_group",groupNR, ".Rdata"))

df <- dfmc; rm(dfmc)

df1 <- subset(df)


df1$value <- (-1)*df1$value
df0 <- df1

input_scores <- unique(df1$score)
plot_folder <- paste0("../Data/Plots/Group ",groupNR,"/")
dir.create(file.path(plot_folder), showWarnings = FALSE)


# drop model == "ens"
df1_save <- df1 
df1 <- subset(df1_save, model != "ens")

# # also drop EMOS.Q and GOF
# df2 <- subset(df1, model != "emos.q" & model != "GOF")
# 
# mypal <- colorspace::rainbow_hcl(10)
# mypal_use <- c("decc.q" = mypal[1],
#                "ecc.q" = mypal[2],
#                "ecc.s" = mypal[3],
#                "gca" = mypal[4],
#                "ssh" = mypal[5],
#                "Clayton" = mypal[6],
#                "Frank" = mypal[7],
#                "Gumbel" = mypal[8],
#                "Surv_Gumbel" = mypal[9],
#                "Seasonal" = mypal[10])
# 
# df2$model <- factor(df2$model, levels = c("decc.q", "ecc.q", "ecc.s", "gca", "ssh", "Clayton", "Frank", "Gumbel","Surv_Gumbel", "Seasonal"))
# model_vec <- c("dECC", "ECC-Q", "ECC-S", "GCA", "SSh", "Clayton", "Frank", "Gumbel","Surv_Gumbel", "Seasonal")

# also drop EMOS.Q, GOF, ECC methods and GCA
df2 <- subset(df1, model != "emos.q" & model != "GOF" & model != "decc.q" & model != "ecc.q" & model != "ecc.s" & model != "gca")

mypal <- colorspace::rainbow_hcl(10)
mypal_use <- c("ssh" = mypal[1],
               "Clayton" = mypal[2],
               "Frank" = mypal[3],
               "Gumbel" = mypal[4],
               "Surv_Gumbel" = mypal[5],
               "Seasonal" = mypal[6])

df2$model <- factor(df2$model, levels = c("ssh", "Clayton", "Frank", "Gumbel","Surv_Gumbel", "Seasonal"))
model_vec <- c("SSh", "Clayton", "Frank", "Gumbel","Surv_Gumbel", "Seasonal")

ylimitFunc <- function(val1, val2) {
  return(1.5 * max(abs(val1), abs(val2)))
}


# Plot the scores for copula

  


plotScores <- function(dfplot, scval) {
  alpha <- 0.25
  
  quants <- unname(quantile(dfplot$value, c(0.01, 0.99)))
  
  ylimits <- c(1.5 * min(quants[1], qnorm(alpha)), 1.5 * max(quants[2], qnorm(1 - alpha)))
  
  
  
  p1 <- ggplot(dfplot, aes(model, value, colour = model))
  p1 <- p1 + ylim(ylimits[1], ylimits[2])
  p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(alpha), ymax=qnorm(1-alpha)), fill = "gray75", color="gray75", alpha=alpha)
  p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  
  
  if(scval == "es_list"){ 
    title <- paste("Energy Score")
    p1 <- p1 + ggtitle(title)
  } 
  if(scval == "vs0_list"){
    title <- paste("Variogram Score")
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

  
plotWidth <- 9
plotHeight <- 6

# Get the plots
p1 <- plotScores(dfplot, this_score)


# Save the plots individually
saveFigure(paste0("ES_group_", groupNR, ".png"), p1)

## VS
this_score <- "vs0_list"
dfplot <- subset(df2, score == this_score)


plotWidth <- 9
plotHeight <- 6

# Get the plots
p1 <- plotScores(dfplot, this_score)


# Save the plots individually
saveFigure(paste0("VS_group_", groupNR, ".png"), p1)

