rm(list=ls())

library(ggplot2)
library(gridExtra)
library(here)

setwd(paste0(here("multiv_pp-master"), "/reproducing results and figures"))

groupNR <- 5
fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable


df <- res; rm(res)

df1 <- df$timing_list

input_times <- names(df1)



dfplot <- data.frame(matrix(ncol=2, nrow=0))

names(dfplot) <- c("input", "value")



for (input_name in input_times) {
  val <- df$timing_list[[input_name]]
  newdf <- data.frame(input = input_name, value = val)
  dfplot <- rbind(dfplot, newdf)
}




plot_folder <- paste0("../Data/Plots/Group ",groupNR,"/")
dir.create(file.path(plot_folder), showWarnings = FALSE)


dfplot <- subset(dfplot, input != "emos.q" & input != "ens" & input != "GOF")


mypal <- colorspace::rainbow_hcl(13)
mypal_use <- c("uvpp" = mypal[3],
               "decc.q" = mypal[5],
               "ecc.q" = mypal[6],
               "ecc.s" = mypal[7],
               "gca" = mypal[8],
               "ssh" = mypal[9],
               "Clayton" = mypal[10],
               "Frank" = mypal[11],
               "Gumbel" = mypal[12],
               "Surv_Gumbel" = mypal[12],
               "Seasonal" = mypal[13])



time_obs <- c("UVPP", "dECC", "ECC-Q", "ECC-S", "GCA", "SSh", "Clayton", "Frank", "Gumbel", "Surv_Gumbel", "Seasonal")

dfplot$input <- factor(dfplot$input, levels = 
                         c("uvpp", "decc.q", "ecc.q", "ecc.s", "gca", "ssh", "Clayton", "Frank", "Gumbel","Surv_Gumbel", "Seasonal"))


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
fileName <- paste0("Times_group_", groupNR, ".png")

ggsave(
  paste0(plot_folder, fileName),
  p1,
  width = plotWidth,
  height = plotHeight,
  dpi = res
)



