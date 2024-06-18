# rm(list=ls())

library(ggplot2)
library(gridExtra)

plot_dm_scores <- function(timeWindow, fix_training_days, training_days_method) {
  for (groupNR in 1:5) {
    fName <- paste0("Res_group_", groupNR)

    if (fix_training_days) {
      fName <- paste0(training_days_method, "_", fName)
    }

    fName <- paste0("m_", timeWindow, "_", fName)

    print(fName)

    load(paste0("Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable

    savedir <- paste0("Data/TestStatistic/")
    savename <- paste0("TestStatistic_data_group", groupNR, ".Rdata")

    if (fix_training_days) {
      savename <- paste0(training_days_method, "_", savename)
    }

    savename <- paste0("m_", timeWindow, "_", savename)

    load(paste0(savedir, savename))

    # Load and transform data
    df1 <- subset(dfmc)
    df1$value <- (-1) * df1$value
    df0 <- df1

    # Extract possible scoring options
    input_scores <- unique(df1$score)
    plot_folder <- paste0("Data/Plots/Group ", groupNR, "/")

    if (fix_training_days) {
      plot_folder <- paste0(plot_folder, training_days_method, "/")
    }

    dir.create(file.path(plot_folder), showWarnings = FALSE)


    # drop model == "ens"
    df1_save <- df1
    df1 <- subset(df1_save, model != "ens")


    # also drop EMOS.Q, GOF, ECC methods and GCA
    df2 <- subset(df1, model != "emos.q" & model != "GOF" & model != "decc.q" & model != "ecc.q" & model != "ecc.s" & model != "ssh")

    mypal <- colorspace::rainbow_hcl(20)
    mypal_use <- c(
      "ssh.h"       	  	= mypal[1],
      "ssh.i"   	  	  	= mypal[2],
      "gca" 	  	  	  	= mypal[3],
      "gca.sh" 	  		  	= mypal[4],
      "gca.cop" 	  	  	= mypal[5],
      "gca.cop.sh"  	  	= mypal[6],
      "Clayton" 	  	  	= mypal[7],
      "Claytonsh" 		  	= mypal[8],
      "Frank" 	  		  	= mypal[9],
      "Franksh" 	  	  	= mypal[10],
      "Gumbel" 	  		  	= mypal[11],
      "Gumbelsh" 	 		  	= mypal[12],
      "Surv_Gumbel" 	  	= mypal[13],
      "Surv_Gumbelsh" 	  = mypal[14]
    )


    # Map technical model name and human friendly model name
    levels <- c("ssh.h", "ssh.i", "gca", "gca.sh", "gca.cop", "gca.cop.sh",
                "Clayton", "Claytonsh", "Frank", "Franksh", "Gumbel", "Gumbelsh", "Surv_Gumbel", "Surv_Gumbelsh")
    model_vec <- c("SSh-H", "SSh-I14", "GCA", "GCAsh", "CopGCA", "CopGCAsh", 
                   "Clayton", "ClaytonSh", "Frank", "FrankSh", "Gumbel", "GumbelSh", "Surv_Gumbel", "Surv_GumbelSh")
    model2display_names <- data.frame(model_names = levels, display_names = model_vec)

    # Map technical scoring name and human friendly scoring name
    scores <- c("es_list", "vs0_list", "crps_list", "crps_1", "crps_2", "crps_3")
    scores_vec <- c("ES", "VS", "CRPS (Mult)", "CRPS (1)", "CRPS (2)", "CRPS (3)")
    score2display_score <- data.frame(score_names = scores, display_names = scores_vec)

    # Change representation with model levels
    df2$model <- factor(df2$model, levels = levels)

    ylimitFunc <- function(val1, val2) {
      return(1.5 * max(abs(val1), abs(val2)))
    }


    # Plot the scores for copula
    plotScores <- function(dfplot, scval) {
      alpha <- 0.25

      quants <- unname(quantile(dfplot$value, c(0.01, 0.99)))

      ylimits <- c(1.5 * min(quants[1], qnorm(alpha)), 1.5 * max(quants[2], qnorm(1 - alpha)))

      # Create the plot
      p1 <- ggplot(dfplot, aes(model, value, colour = model))
      p1 <- p1 + ylim(ylimits[1], ylimits[2])
      p1 <- p1 + geom_rect(mapping = aes(xmin = -Inf, xmax = Inf, ymin = qnorm(alpha), ymax = qnorm(1 - alpha)), fill = "gray75", color = "gray75", alpha = alpha)
      p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")

      # Styling
      title <- paste(subset(score2display_score, score_names == score)$display_names)
      p1 <- p1 + ggtitle(title)

      # X labels
      model_vec <- c()
      for (model in sort(dfplot$model)) {
        model_vec <- c(model_vec, subset(model2display_names, model_names == model)$display_names)
      }

      p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
      p1 <- p1 + xlab("Model") + ylab("DM test statistic")
      p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec)
      p1 <- p1 + scale_x_discrete(label = model_vec)

      return(p1)
    }

    # Function to save a plot
    saveFigure <- function(fileName, fig, plotWidth, plotHeight) {
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



    createAndSave <- function(this_score) {
      # Plot data
      dfplot <- subset(df2, score == this_score)
      
      .GlobalEnv$dfplot <- dfplot

      # Plot size
      plotWidth <- 12
      plotHeight <- 6

      # Get the plots
      p1 <- plotScores(dfplot, this_score)


      # Save the plots individually
      saveName <- paste0(subset(score2display_score, score_names == score)$display_names, "_group_", groupNR, ".png")

      if (timeWindow != 0) {
        saveName <- paste0("m_", timeWindow, "_", saveName)
      }

      # Save the plot
      saveFigure(saveName, p1, plotWidth, plotHeight)
    }

    ## Save and plot each score
    for (score in scores)
    {
      createAndSave(score)
    }
  }
}

# for (fix_training_days in c(TRUE, FALSE)) {
#   if (fix_training_days) {
#     for (training_days_method in c("random_past", "last_m_days", "random_2w_interval")) {
#       plot_dm_scores(100, fix_training_days, training_days_method)
#     }
#   } else {
#     plot_dm_scores(100, fix_training_days, training_days_method)
#   }
# }

plot_dm_scores(50, FALSE, training_days_method)
