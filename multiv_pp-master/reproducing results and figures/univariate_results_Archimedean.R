rm(list=ls())

library(ggplot2)
library(gridExtra)

setwd("C:/Users/20192042/OneDrive - TU Eindhoven/Courses/BEP - BAM/Code/multiv_pp-master/reproducing results and figures")

load("../Data/TestStatistic/testStatistic_Archimedean_1.Rdata")

df <- dfmc; rm(dfmc)
df1 <- subset(df, model == "ens" | model == "ecc.q" | model == "gca")
df1 <- subset(df1, score == "crps_list")
df1 <- subset(df1, eps == 1)

df1$value = (-1)*df1$value

df1$model <- factor(df1$model, levels = c("ens", "decc.q", "ecc.q", "ecc.s", "gca", "ssh"))

mypal <- c(colorspace::rainbow_hcl(5), "black")
mypal_use <- c("ens" = mypal[6],
               "ecc.q" = mypal[2],
               "gca" = mypal[4])
model_vec <- c("Ens.", "ECC-Q", "GCA")

p1 <- ggplot(df1, aes(model, value, colour = model)) 
p1 <- p1 + scale_x_discrete(label = model_vec) # label = model_vec
p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
p1 <- p1 + geom_rect(mapping=aes(xmin=-Inf, xmax=Inf, ymin=qnorm(0.025), ymax=qnorm(0.975)), fill = "gray75", color="gray75")
p1 <- p1 + geom_hline(yintercept = 0, linetype = "dashed") 
p1 <- p1 + ggtitle("Setting 1, " * epsilon ~ "= 1") + geom_boxplot(outlier.shape = NA)
p1 <- p1 + xlab("Model") + ylab("DM test statistic") 
p1


g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(p1)

p1_save <- p1 + theme_bw() + theme(legend.position = "none")

p1_save

png("../Data/Plots/crps_Archimedean.png", width = 15, height = 5, units = "in", res = 400, pointsize = 11)
grid.arrange(p1_save, legend,
             ncol = 1,
             widths = c(3))

dev.off()
