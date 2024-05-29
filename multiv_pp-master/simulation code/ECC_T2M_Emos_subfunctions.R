
## CODE IS ADAPTED FROM ##

# -------------------------------------------------------------------
# - NAME:        ECC_T2M_Emos_subfunctions.R
# - AUTHOR:      Perrone Elisa
# - DATE:        2021-02-18
# -------------------------------------------------------------------
# - DESCRIPTION: EMOS -- optimization based on CRPS
# -------------------------------------------------------------------
# - CONTENT:    -> emos_T2M_mean_singleForecast             
#                  (tailored for ALADIN-LAEF data)
# -------------------------------------------------------------------


#####################################################################
#  Function that applies EMOS to a dataset data.lt
#  leadtime and stations are fixed 
#  (adjusted from Moritz's code)
#####################################################################
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


setwd(here2())
source("emos_T2M_mean_singleForecast_subfunctions-orig.R")

emos_T2M_mean_singleForecast<-function(data.lt, trainingDays){
   
   data.lt$ensmean <-   apply(data.lt[,grep('^laef',names(data.lt))],1,mean)
   data.lt$enssd <-     apply(data.lt[,grep('^laef',names(data.lt))],1,sd)
   data.lt$ensvar <-    apply(data.lt[,grep('^laef',names(data.lt))],1,var)
   data.lt$ensmedian <- apply(data.lt[,grep('^laef',names(data.lt))],1,median)

   nd <- length(unique(data.lt$td))
   # dim(data.lt) = c(nd, 23)
   # ALADIN-LAEF has 17 ensemble members, which are in data.lt[, 7:23]
   # the observations are stored in data.lt[,6]
   
   # Could be 50 --> pick same for mvpp
   td.emos<- trainingDays #training days
   npred <-1

   # model estimation
   slopeCoef<- array(dim=c((nd-td.emos), (npred+1)))
   varCoef<-   array(dim=c((nd-td.emos), 2))
   valid_df <- data.frame(matrix(nrow=(nd-td.emos),ncol= 6))

   #loop over days
   for (i.d in 1:(nd-td.emos)){
     
      train.index <- data.lt$ytime %in% data.lt$ytime[i.d:(i.d+td.emos-1)]
      valid.index <- data.lt$ytime %in% data.lt$ytime[(i.d+td.emos)]

      
      # fit linear regression for initial values
      lm.fit <- lm(obs~ensmean, data=data.lt[train.index,])

      #(3b) use initial values according to R Package of Gneiting et al. 2005
      a.init <- coef(lm.fit)[1]
      B.init <- coef(lm.fit)[-1]
      B.init <- sqrt(abs(B.init))
      c.init <- 5
      c.init <- sqrt(c.init)
      d.init <- 1
      d.init <- sqrt(d.init)

      par_init <- c(a.init,B.init,c.init,d.init)

      t.obs <- data.lt$obs[train.index]
      t.ens <- cbind(1,as.matrix(data.lt$ensmean[train.index]))
      t.var <- data.lt$ensvar[train.index]


      #(3c) fitting of parameters using CRPS minimization
      fits <- optim(par_init, crps.normal,
                   t.obs = t.obs,
                   t.ens = t.ens,
                   t.var = t.var,
                   K = dim(t.ens)[2],
                   method = "BFGS",
                   control   = list(maxit=as.integer(1e7))
                   )

      #  (3d) transform parameters to observable scale and store them in data.frame
      fits$par[-1] <-    fits$par[-1]^2
      slopeCoef[i.d,] <- fits$par[1:(npred+1)]
      varCoef[i.d,] <-   fits$par[(npred+2):length(fits$par)]

      #  (4) VALIDATION
      #  (4a) calculate mean, variance, standard deviation
      ens_mu <-  slopeCoef[i.d,]%*%c(1,data.lt$ensmean[valid.index])
      ens_var <- varCoef[i.d,]%*%c(1,data.lt$ensvar[valid.index])
      ens_sd <-  sqrt(ens_var)

      #  (4b) calculate crps
      crps.emos <-    crps.norm(data.lt$obs[valid.index],ens_mu,ens_sd) 
      crps.raw_ens <- crps.ensemble(as.numeric(data.lt$obs[valid.index]),as.matrix(data.lt[valid.index,grep('^laef',names(data.lt))]))

      #  (4c)  store scores in data.frame
      valid_df[i.d,] <-  cbind(data.lt$td[valid.index], data.lt$obs[valid.index], ens_mu, ens_sd, crps.emos, crps.raw_ens)
      names(valid_df) <- c("td", "obs", "ens_mu", "ens_sd", "crps_emos", "crps.raw_ens")

   }# end of loop over days

   return(valid_df)

}

