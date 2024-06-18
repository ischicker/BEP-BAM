#rm(list=ls())

library(ggplot2)

plot_boxplots <-function(timeWindow, fix_training_days, training_days_method) {
  
  for (groupNR in 1:5) {
    
    fName <- paste0("Res_group_", groupNR)
    
    if (fix_training_days) {
      fName <- paste0(training_days_method, "_", fName)
    }
    
    if (timeWindow != 0) {
      fName <- paste0("m_", timeWindow, "_", fName)
    }
    
    load(paste0("Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable
    
    print(fName)
    
    
    
    plot_folder <- paste0("Data/Plots/Group ",groupNR,"/")
    
    if (fix_training_days) {
      plot_folder <- paste0(plot_folder, training_days_method, "/")
    }
    
    dir.create(file.path(plot_folder), showWarnings = FALSE)
    

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
    
    
    input_models <- levels
    input_scores <- c("es_list","vs1_list","vs1w_list","vs0_list","vs0w_list", "crps_list", "crps_1", "crps_2", "crps_3")
    
    
    for (this_score in input_scores) {
      
      dfplot <- data.frame(matrix(ncol=2, nrow=0))
      
      names(dfplot) <- c("input", "value")
      
      
      
      for (input_name in input_models) {
        if(this_score %in% c("crps_1", "crps_2", "crps_3")){
          n <- as.numeric(gsub("[^0-9]", "", this_score))
          val <- c(res[["crps_list"]][[input_name]][,n])
        } else {
          val <- c(res[[this_score]][[input_name]])
        }
        newdf <- data.frame(input = input_name, value = val)
        dfplot <- rbind(dfplot, newdf)
      }
      
      dfplot$input <- factor(dfplot$input, levels = input_models)
      
      
      p1 <- ggplot(dfplot, aes(x = input, y = value, colour = input))
      p1 <- p1 + geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed", color = "gray25")
      
      
      p1 <- p1 + ggtitle(this_score) 
      
      model_vec <- c()
      for (model in input_models) {
        model_vec <- c(model_vec, subset(model2display_names, model_names == model)$display_names)
      }
      
      p1 <- p1 + theme_bw() + theme(legend.position = "bottom")
      p1 <- p1 + xlab("Model") + ylab("Score value") 
      p1 <- p1 + scale_color_manual(values = mypal_use, name = "Model", label = model_vec) 
      p1 <- p1 + scale_x_discrete(label = model_vec)+
        theme(plot.title = element_text(hjust = 0.5))
      
      
      plotWidth <- 12
      plotHeight <- 6
      resolution <- 400
      fileName <- paste0("scores_", this_score, "_group_", groupNR, ".png")
      
      if (timeWindow != 0) {
        fileName <- paste0("m_", timeWindow, "_", fileName)
      }
      
      ggsave(
        paste0(plot_folder, fileName),
        p1,
        width = plotWidth,
        height = plotHeight,
        dpi = resolution
      )
      
    }
  }
}


# for (fix_training_days in c(TRUE, FALSE)) {
#   if (fix_training_days) {
#     for (training_days_method in c("random_past", "last_m_days", "random_2w_interval")) {
#       
#       plot_boxplots(100, fix_training_days, training_days_method)
#       
#     }
#   } else {
#     plot_boxplots(100, fix_training_days, training_days_method)
#   }
# }

plot_boxplots(50, FALSE, training_days_method)
