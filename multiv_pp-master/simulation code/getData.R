library("here")

getData <- function() {
  
  ### Input consists of
  # 72 different lead times
  # 216 stations (input$leadtimen consists of the stations)
  # 17 ensemble members (laef1, ..., leaf17)
  # input$leadtimen$lt contains the lead time n
  # 990 different time stamps (unique values of "initial")
  # td contains the discretized time values
  # initial is the time of measurement
  # ytime = initial + 1
  # inca (Integrated Nowcasting through Comprehensive Analyses) contains the observations
  # cosmo (Consortium for Small-scale Modeling) another forecast

  
  # Data file names
  fName2013 <- "INPUT-DATA_temp_2013071000-2016033000"
  fName2016 <- "INPUT-DATA_temp_2016040100-2018050100"
  
  # Load data separately
  load(paste0("../../Data/", fName2013, ".Rdata")) # loads data in "input" variable
  dat2013 <- input; rm(input)
  
  load(paste0("../../Data/", fName2016, ".Rdata")) # loads data in "input" variable
  dat2016 <- input; rm(input)
  
  # Focus on one specific lead time
  leadTime <- 24
  data2013 <- dat2013[[paste0("leadtime", leadTime)]]
  data2016 <- dat2016[[paste0("leadtime", leadTime)]]
  
  # Delete unnecessary data
  rm(dat2013, dat2016)
  
  # Data from 2013 - 2016 contains cosmo column, whereas 2016 - 2018 does not
  data2013$cosmo <- NULL
  
  # Recreate discrete time steps
  tdmax <- max(data2013$td)
  data2016$td <- data2016$td + tdmax
  
  # Merge data files and remove unused structures
  data_unfiltered <- rbind(data2013, data2016)
  rm(data2013, data2016)
  
  # First part of data contains incorrect LAEF values --> find index
  i <- 1
  while (i <= dim(data_unfiltered)[1]) {
    if (data_unfiltered[i,]$laef1 != 273.15) {
      break
    }
    i <- i + 1
  }
  
  data <- data_unfiltered[i:dim(data_unfiltered)[1],]
  rm(data_unfiltered)
  
  # Time format is YYYYMMDD00 and is converted to Date object
  data$time <- as.Date(sapply(data$initial, FUN = function(x) substr(toString(x), start = 1, stop = 8)), "%Y%m%d")
  
  # Remove entries with NA values and rename one column
  data$obs <- data$inca
  data$inca <- NULL
  data <- data[!unname(apply(data, MARGIN = 1, FUN = function(x) any(is.na(x)))),]
  
  #vector of all distinct stations
  statNR <-  unique(data$stat) 
  
  # Groups in Perrone et al. 2020 are:
  # group 1: stations 188 - 189 - 191
  # group 2: stations 64 - 169 - 188
  # group 3: stations 19 - 150 - 172
  group1 <- c(188, 189, 191)
  group2 <- c(64, 169, 188)
  group3 <- c(19, 150, 172)
  group4 <- match(c(11024, 11077, 11022), statNR)
  group5 <- match(c(11007, 11036, 11056), statNR)
  
  # Select the data for those groups
  data1 <- subset(data, stat %in% statNR[group1])
  data2 <- subset(data, stat %in% statNR[group2])
  data3 <- subset(data, stat %in% statNR[group3])
  data4 <- subset(data, stat %in% statNR[group4])
  data5 <- subset(data, stat %in% statNR[group5])
  
  
  # Some time stamps do not have data for all stations (after removing NA) --> remove data points
  cleanup <- function(dat){
    tds <- unique(dat$td)
    stats <- length(unique(dat$stat))
    
    # print(tds)
    
    for (td.val in tds) {
      # Imbalance
      if (sum(dat$td == td.val) != stats) {
        dat <- subset(dat, td != td.val)
      }
    }
    return(dat)
  }
  
  
  data1 <- cleanup(data1)
  data2 <- cleanup(data2)
  data3 <- cleanup(data3)
  data4 <- cleanup(data4)
  data5 <- cleanup(data5)
  
  .GlobalEnv$statNR <- statNR
  
  .GlobalEnv$group1 <- group1
  .GlobalEnv$group2 <- group2
  .GlobalEnv$group3 <- group3
  .GlobalEnv$group4 <- group4
  .GlobalEnv$group5 <- group5
  
  .GlobalEnv$data1 <- data1
  .GlobalEnv$data2 <- data2
  .GlobalEnv$data3 <- data3
  .GlobalEnv$data4 <- data4
  .GlobalEnv$data5 <- data5
  
  .GlobalEnv$data <- data
}