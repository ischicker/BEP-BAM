# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(zoo)

here2 <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  if ("RStudio" %in% args) {
    dirname(rstudioapi::getActiveDocumentContext()$path)
  } else {
    file_arg <- "--file="
    filepath <- sub(file_arg, "", grep(file_arg, args, value = TRUE))
    dirname(filepath)
  }
}

# Function to create and save plots for a single variable with optional rolling average
create_plots <- function(data, var_name, rolling_avg = FALSE, window_size = 10, save_dir = "plots") {
  
  if (!is.null(dim(data)) && dim(data)[1] == 1) {
    data <- data[1,]
  } 
  
  
  
  
  if (!is.null(dim(data))) {
    
    if (rolling_avg) {
      data <- apply(data, 2, function(x) rollmean(x, window_size, fill = NA))
    }
    
    rows <- dim(data)[1]
    columns <- dim(data)[2]
  
    # Create a dataframe for ggplot
    df <- data.frame(
      Index = rep(1:rows, 3),
      Value = c(data[, 1], data[, 2], data[, 3]),
      Column = factor(rep(1:3, each = rows))
    )
    
    # Generate plots for each column
    p1 <- ggplot(df[df$Column == 1, ], aes(x = Index, y = Value)) +
      geom_line() +
      ggtitle(paste(var_name, " - Station 1"))
    
    p2 <- ggplot(df[df$Column == 2, ], aes(x = Index, y = Value)) +
      geom_line() +
      ggtitle(paste(var_name, " - Station 2"))
    
    p3 <- ggplot(df[df$Column == 3, ], aes(x = Index, y = Value)) +
      geom_line() +
      ggtitle(paste(var_name, " - Station 3"))
    
    # Combine the plots into one
    combined_plot <- grid.arrange(p1, p2, p3, ncol = 1)
  } else {
    
    
    
    # Convert the data to a dataframe
    df <- as.data.frame(data)
    colnames(df) <- c("y")
    
    if (rolling_avg) {
      df <- apply(df, 2, function(x) rollmean(x, window_size, fill = NA))
    }
    
    # Generate the plot
    combined_plot <- ggplot(df, aes(x = 1:nrow(df), y = y)) +
      geom_line() +
      ggtitle(var_name) +
      xlab("Index") +
      ylab("Value")
    
  }
  
  # Create the directory if it does not exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir)
  }
  
  if (rolling_avg) {
    savename <- paste0(var_name, "_rm_", window_size, "_scores.png")
  } else {
    savename <- paste0(var_name, "_scores.png")
  }
  
  # Save the combined plot
  ggsave(filename = file.path(save_dir, savename),
         plot = combined_plot,
         width = 16, height = 12, dpi = 300)
}

setwd(here2())


groupNR <- 3
fName <- paste0("Res_group_", groupNR)
load(paste0("../Data/Rdata_LAEF/", fName, ".Rdata")) # loads data in "res" variable

scores <- names(res)[1:6]

for (rolling_avg in c(TRUE, FALSE)) {
  for (score in scores) {
    # Specify the directory to save the plots
    save_directory <- paste0("../Data/Plots/Group ",groupNR,"/", sub("_list$", "", score), "/")
    
    # Loop through each variable in the list and create plots
    data <- res[[score]]
    for (var_name in names(data)) {
      create_plots(data[[var_name]], var_name, rolling_avg = rolling_avg, window_size = 30, save_dir = save_directory)
    }
  }
}


