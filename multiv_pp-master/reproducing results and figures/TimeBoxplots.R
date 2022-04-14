rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")

# Model 1 : Standard Gaussian Marginals

modelSetting <- 1
setting <- 2

fName <- paste0("Archimedean","_setting_",setting,"_model_", modelSetting,"_ID_")


flist <- list.files("../Data/Rdata/")
existing <- as.numeric(sapply(flist, FUN = function(x) as.numeric(strsplit(strsplit(x, fName)[[1]][2], ".Rdata"))))


load(paste0("../Data/Rdata/",fName, "1.Rdata"))

df <- res; rm(res)

df1 <- df$timing_list

input_times <- names(df1)


dfplot <- data.frame(matrix(ncol=2, nrow=0))

names(dfplot) <- c("input", "value")


for(ID in existing[!is.na(existing)]){
  load(paste0("../Data/Rdata/", fName, ID,".Rdata"))
  avg_times <- lapply(res$timing_list, mean )
  for (input_name in input_times) {
    val <- avg_times[[input_name]]
    newdf <- data.frame(input = input_name, value = val)
    dfplot <- rbind(dfplot, newdf)
  }
  
}





plot_folder <- paste0("../Data/Plots/Arch_setting_",setting,"_model_",modelSetting,"/")
dir.create(file.path(plot_folder), showWarnings = FALSE)


dfplot <- subset(dfplot, input != "emos.q")


mypal <- colorspace::rainbow_hcl(12)
mypal_use <- c("obs" = mypal[1],
               "fc" = mypal[2],
               "uvpp" = mypal[3],
               "ens" = mypal[4],
               "decc.q" = mypal[5],
               "ecc.q" = mypal[6],
               "ecc.s" = mypal[7],
               "gca" = mypal[8],
               "ssh" = mypal[9],
               "clayton" = mypal[10],
               "frank" = mypal[11],
               "gumbel" = mypal[12])



time_obs <- c("Observations", "Forecast", "UVPP", "Ens", "dECC", "ECC-Q", "ECC-S", "GCA", "SSh", "Clayton", "Frank", "Gumbel")

dfplot$input <- factor(dfplot$input, levels = 
                         c("obs", "fc", "uvpp", "ens", "decc.q", "ecc.q", "ecc.s", "gca", "ssh", "clayton", "frank", "gumbel"))


p1 <- ggplot(dfplot, aes(x = input, y = value, colour = input))
p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")


p1 <- p1 + ggtitle(paste0("Timing Boxplots")) 



p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
p1 <- p1 + xlab("Action") + ylab("Time (s)") 
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Action", label = time_obs) 
p1 <- p1 + scale_x_discrete(label = time_obs)+
  theme(plot.title = element_text(hjust = 0.5))


plotWidth <- 9
plotHeight <- 6
res <- 400
fileName <- paste0("Arch_Times_setting_", setting, "_model_",modelSetting,".png")

ggsave(
  paste0(plot_folder, fileName),
  p1,
  width = plotWidth,
  height = plotHeight,
  dpi = res
)



