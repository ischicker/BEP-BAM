rm(list=ls())

library(ggplot2)
library("here")


setwd(paste0(here("multiv_pp-master"), "/reproducing results and figures"))


groupNR <- 1

fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable


input_models <- c("ssh","gca","Clayton","Frank","Gumbel","Surv_Gumbel","Seasonal")
input_scores <- c("es_list","vs1_list","vs1w_list","vs0_list","vs0w_list")


plot_folder <- paste0("../Data/Plots/Group ",groupNR,"/")
dir.create(file.path(plot_folder), showWarnings = FALSE)




mypal <- colorspace::rainbow_hcl(7)
mypal_use <- c("ssh" = mypal[1],
               "gca" = mypal[2],
               "Clayton" = mypal[3],
               "Frank" = mypal[4],
               "Gumbel" = mypal[5],
               "Surv_Gumbel" = mypal[6],
               "Seasonal" = mypal[7])



x_labels <- c("SSh", "GCA", "Clayton", "Frank", "Gumbel", "Surv_Gumbel", "Seasonal")


for (this_score in input_scores) {
  
  dfplot <- data.frame(matrix(ncol=2, nrow=0))
  
  names(dfplot) <- c("input", "value")
  
  
  
  for (input_name in input_models) {
    val <- c(res[[this_score]][[input_name]])
    newdf <- data.frame(input = input_name, value = val)
    dfplot <- rbind(dfplot, newdf)
  }
  
  dfplot$input <- factor(dfplot$input, levels = input_models)
  
  
  p1 <- ggplot(dfplot, aes(x = input, y = value, colour = input))
  p1 <- p1 + geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
  
  
  p1 <- p1 + ggtitle(this_score) 
  
  
  
  p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
  p1 <- p1 + xlab("Model") + ylab("Score value") 
  p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = x_labels) 
  p1 <- p1 + scale_x_discrete(label = x_labels)+
    theme(plot.title = element_text(hjust = 0.5))
  
  
  plotWidth <- 9
  plotHeight <- 6
  resolution <- 400
  fileName <- paste0("Scores_", this_score, "_group_", groupNR, ".png")
  
  ggsave(
    paste0(plot_folder, fileName),
    p1,
    width = plotWidth,
    height = plotHeight,
    dpi = resolution
  )
  
}

