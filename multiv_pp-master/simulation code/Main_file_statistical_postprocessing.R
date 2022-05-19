##################################################
#### Parametric Rscript for temperature  #########
####       Elisa Perrone  2018           #########
##################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=7) {
  stop("Usage: Rscript main_3st.R [output_name] [current_dir] [leadtime] [runs] [station1] [station2] [station3]", call.=FALSE)
}

print(Sys.time())

library("stats")
library("compiler")
library('rlist')

add_name <- args[1]
main_path <- args[2]
ld.t<-as.double(args[3])
n.runs<-as.double(args[4])
stations<-c(as.double(args[5]),as.double(args[6]),as.double(args[7]))

k1<-stations[1]
k2<-stations[2]
k3<-stations[3]

# The groups I considered in Perrone et al. 2020 are:
# group 1: stations 188 - 189 - 191
# group 2: stations 64 - 169 - 188
# group 3: stations 19 - 150 - 172

npred <-1     
td.emos<-50
sample.size<-17


input_path <- main_path
output_path <-  paste0(main_path, "/outcome_files/")
output_file <-  paste0(output_path,"/",add_name,"_ldt",as.character(ld.t),"_Runs",as.character(n.runs),"_StNr",as.character(k1),"_",as.character(k2),"_",as.character(k3),".Rdata")
data_path<-"/Users/elisa_perrone/Dropbox/My Mac (C02YRA8ZLVCG)/Desktop/ZAMG/Stuff_Internship/ZAMG_saved/"

source(paste0(input_path, "/emos_T2M_mean_singleForecast_subfunctions-orig.R"))
source(paste0(input_path, "/ECC_T2M_emos_subfunctions copy.R"))



# -------------------------------------------------------------------
#  (2) DATA READING
# -------------------------------------------------------------------
load(paste0(data_path,"/INPUT-DATA_temp_2013071000-2016033000.Rdata"))
#select 30 months of data
input_test1<- lapply(input, FUN = function(x) subset(x,initial>=2014020100 & initial<=2016033000))
input_test1<-input_test1[ld.t] #we select the leadtime 12
statNR <-  input$leadtime1$stat[1:216] #vector of stations


load(paste0(data_path,"/INPUT-DATA_temp_2016040100-2018050100.Rdata"))

if(ld.t<24){
  input_test2<- lapply(input, FUN = function(x) subset(x,ytime>=2016040100 & ytime<=2018050100))
}else{
  input_test2<- lapply(input, FUN = function(x) subset(x,initial>=2016040100 & initial<=2018050100))
}
input_test2<-input_test2[ld.t] #we select the leadtime 12
statNR <-  input$leadtime1$stat[1:216] #vector of stations
rm(input)#remove unnecessary data

input_test<-list(rbind(input_test1[[1]][1:23],input_test2[[1]]))


# -------------------------------------------------------------------
#  (2b) Selection of past observations
# -------------------------------------------------------------------

#select the records according to statNR and the chosen leadtime 12
data.stat <- lapply(input_test, FUN = function(x) subset(x,stat==statNR[k1]))
data.ltN1 <-   data.stat[[1]]

#prepare data (calculate means, sd, ...)
names(data.ltN1)[(names(data.ltN1)=="inca")] <- "obs"
stN1<-na.omit(data.ltN1)
#stN1<-stN1[!duplicated(stN1$obs),]

#select the records according to statNR and the chosen leadtime 12
data.stat <- lapply(input_test, FUN = function(x) subset(x,stat==statNR[k2]))
data.ltN2 <-   data.stat[[1]]

#prepare data (calculate means, sd, ...)
names(data.ltN2)[(names(data.ltN2)=="inca")] <- "obs"
stN2<-na.omit(data.ltN2)
#stN2<-stN2[!duplicated(stN2$obs),]

#select the records according to statNR and the chosen leadtime 12
data.stat <- lapply(input_test, FUN = function(x) subset(x,stat==statNR[k3]))
data.ltN3 <-   data.stat[[1]]

#prepare data (calculate means, sd, ...)
names(data.ltN3)[(names(data.ltN3)=="inca")] <- "obs"
stN3<-na.omit(data.ltN3)
#stN3<-stN3[!duplicated(stN3$obs),]

length(stN1$ytime)

stN1_N2<-intersect(stN1$ytime,stN2$ytime)
stN2_N3<-intersect(stN2$ytime,stN3$ytime)
stN1_N2_N3<-intersect(stN1_N2,stN2_N3)

length(stN1_N2_N3)

r1<-array("FALSE",dim=length(data.ltN1[,3]))

for(j in 1:length(stN1_N2_N3)){
  for(i in 1:length(data.ltN1[,3])){
    if(stN1_N2_N3[j]==data.ltN1[i,3]){
      r1[i]<-"TRUE"
    }
  }
}
f1<-which(r1 =="FALSE")
data.ltN1<-data.ltN1[-f1,]
dim(data.ltN1)

r2<-array("FALSE",dim=length(data.ltN2[,3]))

for(j in 1:length(stN1_N2_N3)){
  for(i in 1:length(data.ltN2[,3])){
    if(stN1_N2_N3[j]==data.ltN2[i,3]){
      r2[i]<-"TRUE"
    }
  }
}
f2<-which(r2 =="FALSE")
data.ltN2<-data.ltN2[-f2,]
dim(data.ltN2)


r3<-array("FALSE",dim=length(data.ltN3[,3]))

for(j in 1:length(stN1_N2_N3)){
  for(i in 1:length(data.ltN3[,3])){
    if(stN1_N2_N3[j]==data.ltN3[i,3]){
      r3[i]<-"TRUE"
    }
  }
}
f3<-which(r3 =="FALSE")
data.ltN3<-data.ltN3[-f3,]
dim(data.ltN3)

M<-sample.size

# -------------------------------------------------------------------
#  (3) Correct the forecasts and get a new sample
# -------------------------------------------------------------------
mean_varN1<-vector("list")
mean_varN2<-vector("list")
mean_varN3<-vector("list")
Mean_Var_past3<-vector("list")

data.ltN1$ensmean <-   apply(data.ltN1[,grep('^laef',names(data.ltN1))],1,mean)
data.ltN1$enssd <-     apply(data.ltN1[,grep('^laef',names(data.ltN1))],1,sd)
data.ltN1$ensvar <-    apply(data.ltN1[,grep('^laef',names(data.ltN1))],1,var)
data.ltN1$ensmedian <- apply(data.ltN1[,grep('^laef',names(data.ltN1))],1,median)

data.ltN2$ensmean <-   apply(data.ltN2[,grep('^laef',names(data.ltN2))],1,mean)
data.ltN2$enssd <-     apply(data.ltN2[,grep('^laef',names(data.ltN2))],1,sd)
data.ltN2$ensvar <-    apply(data.ltN2[,grep('^laef',names(data.ltN2))],1,var)
data.ltN2$ensmedian <- apply(data.ltN2[,grep('^laef',names(data.ltN2))],1,median)

data.ltN3$ensmean <-   apply(data.ltN3[,grep('^laef',names(data.ltN3))],1,mean)
data.ltN3$enssd <-     apply(data.ltN3[,grep('^laef',names(data.ltN3))],1,sd)
data.ltN3$ensvar <-    apply(data.ltN3[,grep('^laef',names(data.ltN3))],1,var)
data.ltN3$ensmedian <- apply(data.ltN3[,grep('^laef',names(data.ltN3))],1,median)

mean_varN1<-list(data.ltN1$ytime,data.ltN1$ensmean,data.ltN1$enssd,data.ltN1$ensvar,data.ltN1$ensmedian)
names(mean_varN1) <- c("ytime","ensmean","enssd", "ensvar", "ensmedian")

mean_varN2<-list(data.ltN2$ytime,data.ltN2$ensmean,data.ltN2$enssd,data.ltN2$ensvar,data.ltN2$ensmedian)
names(mean_varN2) <- c("ytime","ensmean","enssd", "ensvar", "ensmedian")

mean_varN3<-list(data.ltN3$ytime,data.ltN3$ensmean,data.ltN3$enssd,data.ltN3$ensvar,data.ltN3$ensmedian)
names(mean_varN3) <- c("ytime","ensmean","enssd", "ensvar", "ensmedian")


Mean_Var_past3<-list(mean_varN1,mean_varN2,mean_varN3)
names(Mean_Var_past3)<-c("St11035","St11216","St11320")

nd<-length(Mean_Var_past3$St11035$ytime)
st3<-c(statNR[k1],statNR[k2],statNR[k3])


# -------------------------------------------------------------------
#  (5) Correct the distribution
# -------------------------------------------------------------------

# apply univariate emos to get ens_mu and ens_sd
valid_dfN1 <- emos_T2M_mean_singleForecast(data.ltN1)
valid_dfN2 <- emos_T2M_mean_singleForecast(data.ltN2)
valid_dfN3 <- emos_T2M_mean_singleForecast(data.ltN3)

# -------------------------------------------------------------------
#  (6) Sample (uniform quantiles)
# -------------------------------------------------------------------

# sample from a normal distribution with parameters ens_mu and ens_sd (uniform quantile)
samQ_norN1<-lapply(1:(nd-td.emos),function(x) qnorm(seq(1/(sample.size+1), (1-1/(sample.size+1)),1/(sample.size+1)),mean=valid_dfN1$ens_mu[x],sd=valid_dfN1$ens_sd[x])) #ECC-Q N1
samQ_norN2<-lapply(1:(nd-td.emos),function(x) qnorm(seq(1/(sample.size+1), (1-1/(sample.size+1)),1/(sample.size+1)),mean=valid_dfN2$ens_mu[x],sd=valid_dfN2$ens_sd[x])) #ECC-Q N2
samQ_norN3<-lapply(1:(nd-td.emos),function(x) qnorm(seq(1/(sample.size+1), (1-1/(sample.size+1)),1/(sample.size+1)),mean=valid_dfN3$ens_mu[x],sd=valid_dfN3$ens_sd[x])) #ECC-Q N3

# save the used sample
names(samQ_norN1) <- sprintf("day%i",1:(nd-td.emos))
names(samQ_norN2) <- sprintf("day%i",1:(nd-td.emos))
names(samQ_norN3) <- sprintf("day%i",1:(nd-td.emos))
emosQ_N1_N2_N3 <-list(as.data.frame(samQ_norN1),as.data.frame(samQ_norN2),as.data.frame(samQ_norN3))      

#sample with obs
samQ_norN1_obs  <-lapply(1:(nd-td.emos), function(x) append(samQ_norN1[[x]], valid_dfN1$obs[x], after = length(samQ_norN1[[x]]))) # EMOS Q
samQ_norN2_obs  <-lapply(1:(nd-td.emos), function(x) append(samQ_norN2[[x]], valid_dfN2$obs[x], after = length(samQ_norN2[[x]]))) # EMOS Q
samQ_norN3_obs  <-lapply(1:(nd-td.emos), function(x) append(samQ_norN3[[x]], valid_dfN3$obs[x], after = length(samQ_norN3[[x]]))) # EMOS Q


# save the used sample + obs
names(samQ_norN1_obs) <- sprintf("day%i",1:(nd-td.emos))
names(samQ_norN2_obs) <- sprintf("day%i",1:(nd-td.emos))
names(samQ_norN3_obs) <- sprintf("day%i",1:(nd-td.emos))
emosQ_N1_N2_N3_obs <-list(as.data.frame(samQ_norN1_obs),as.data.frame(samQ_norN2_obs),as.data.frame(samQ_norN3_obs))


# -------------------------------------------------------------------
#  (7) Sample (n.runs random samples)
# -------------------------------------------------------------------

emos_N1_N2_N3<-list()
emos_N1_N2_N3_obs<-list()

for(i in 1:n.runs){
  
  # sample from a normal distribution with parameters ens_mu and ens_sd (uniform quantile)
  sam_norN1<-lapply(1:(nd-td.emos),function(x) rnorm(sample.size, mean=valid_dfN1$ens_mu[x],sd=valid_dfN1$ens_sd[x])) #ECC-Q N1
  sam_norN2<-lapply(1:(nd-td.emos),function(x) rnorm(sample.size, mean=valid_dfN2$ens_mu[x],sd=valid_dfN2$ens_sd[x])) #ECC-Q N2
  sam_norN3<-lapply(1:(nd-td.emos),function(x) rnorm(sample.size, mean=valid_dfN3$ens_mu[x],sd=valid_dfN3$ens_sd[x])) #ECC-Q N3
  
  # save the used sample
  names(sam_norN1) <- sprintf("day%i",1:(nd-td.emos))
  names(sam_norN2) <- sprintf("day%i",1:(nd-td.emos))
  names(sam_norN3) <- sprintf("day%i",1:(nd-td.emos))
  emos_N1_N2_N3[[i]] <-list(as.data.frame(sam_norN1),as.data.frame(sam_norN2),as.data.frame(sam_norN3))
  
  #sample with obs
  sam_norN1_obs  <-lapply(1:(nd-td.emos), function(x) append(sam_norN1[[x]], valid_dfN1$obs[x], after = length(sam_norN1[[x]]))) # EMOS random
  sam_norN2_obs  <-lapply(1:(nd-td.emos), function(x) append(sam_norN2[[x]], valid_dfN2$obs[x], after = length(sam_norN2[[x]]))) # EMOS random
  sam_norN3_obs  <-lapply(1:(nd-td.emos), function(x) append(sam_norN3[[x]], valid_dfN3$obs[x], after = length(sam_norN3[[x]]))) # EMOS random
  
  
  # save the used sample + obs
  names(sam_norN1_obs) <- sprintf("day%i",1:(nd-td.emos))
  names(sam_norN2_obs) <- sprintf("day%i",1:(nd-td.emos))
  names(sam_norN3_obs) <- sprintf("day%i",1:(nd-td.emos))
  emos_N1_N2_N3_obs[[i]] <-list(as.data.frame(sam_norN1_obs),as.data.frame(sam_norN2_obs),as.data.frame(sam_norN3_obs))

}


# -------------------------------------------------------------------
#  (8) ECC-R step 
# -------------------------------------------------------------------
s<-dim(as.data.frame(sam_norN1_obs))[2]

ECCstN1_N2_N3_obs<-list()
EmosN1_N2_N3_obs<-list()


for(k in 1:n.runs){

ECC_stN1<-array(dim=c(s, sample.size))
ECC_stN2<-array(dim=c(s, sample.size))
ECC_stN3<-array(dim=c(s, sample.size))
#
ECC_stN1_obs<-array(dim=c(s, sample.size+1))
ECC_stN2_obs<-array(dim=c(s, sample.size+1))
ECC_stN3_obs<-array(dim=c(s, sample.size+1))

Emos_stN1<-array(dim=c(s,sample.size))
Emos_stN2<-array(dim=c(s,sample.size))
Emos_stN3<-array(dim=c(s,sample.size))
#
Emos_stN1_obs<-array(dim=c(s,sample.size+1))
Emos_stN2_obs<-array(dim=c(s,sample.size+1))
Emos_stN3_obs<-array(dim=c(s,sample.size+1))


for(i in 1:s){
  #  print(i)
  ECC_stN1[i,]<-order.vec(data.ltN1[td.emos + i,7:23],emos_N1_N2_N3[[k]][[1]][[i]])
  ECC_stN2[i,]<-order.vec(data.ltN2[td.emos + i,7:23],emos_N1_N2_N3[[k]][[2]][[i]])
  ECC_stN3[i,]<-order.vec(data.ltN3[td.emos + i,7:23],emos_N1_N2_N3[[k]][[3]][[i]])
  
  ECC_stN1_obs[i,]<- c(ECC_stN1[i,],emos_N1_N2_N3_obs[[k]][[1]][[i]][18])
  ECC_stN2_obs[i,]<-c(ECC_stN2[i,],emos_N1_N2_N3_obs[[k]][[2]][[i]][18])
  ECC_stN3_obs[i,]<-c(ECC_stN3[i,],emos_N1_N2_N3_obs[[k]][[3]][[i]][18])
  
  Emos_stN1[i,]<-emos_N1_N2_N3[[k]][[1]][[i]]
  Emos_stN2[i,]<-emos_N1_N2_N3[[k]][[2]][[i]]
  Emos_stN3[i,]<-emos_N1_N2_N3[[k]][[3]][[i]]
  
  Emos_stN1_obs[i,]<- c(Emos_stN1[i,],emos_N1_N2_N3_obs[[k]][[1]][[i]][18])
  Emos_stN2_obs[i,]<-c(Emos_stN2[i,],emos_N1_N2_N3_obs[[k]][[2]][[i]][18])
  Emos_stN3_obs[i,]<-c(Emos_stN3[i,],emos_N1_N2_N3_obs[[k]][[3]][[i]][18])
  
}

ECCstN1_N2_N3_obs[[k]]<-list(ECC_stN1_obs,ECC_stN2_obs,ECC_stN3_obs)

EmosN1_N2_N3_obs[[k]]<-list(Emos_stN1_obs,Emos_stN2_obs,Emos_stN3_obs)

}

# -------------------------------------------------------------------
#  (9) ECC-Q step 
# -------------------------------------------------------------------


ECCQ_stN1<-array(dim=c(s, sample.size))
ECCQ_stN2<-array(dim=c(s, sample.size))
ECCQ_stN3<-array(dim=c(s, sample.size))
#
ECCQ_stN1_obs<-array(dim=c(s, sample.size+1))
ECCQ_stN2_obs<-array(dim=c(s, sample.size+1))
ECCQ_stN3_obs<-array(dim=c(s, sample.size+1))

EmosQ_stN1<-array(dim=c(s,sample.size))
EmosQ_stN2<-array(dim=c(s,sample.size))
EmosQ_stN3<-array(dim=c(s,sample.size))
#
EmosQ_stN1_obs<-array(dim=c(s,sample.size+1))
EmosQ_stN2_obs<-array(dim=c(s,sample.size+1))
EmosQ_stN3_obs<-array(dim=c(s,sample.size+1))

for(i in 1:s){
  #  print(i)
  ECCQ_stN1[i,]<-order.vec(data.ltN1[td.emos + i,7:23],emosQ_N1_N2_N3[[1]][[i]])
  ECCQ_stN2[i,]<-order.vec(data.ltN2[td.emos + i,7:23],emosQ_N1_N2_N3[[2]][[i]])
  ECCQ_stN3[i,]<-order.vec(data.ltN3[td.emos + i,7:23],emosQ_N1_N2_N3[[3]][[i]])
  
  ECCQ_stN1_obs[i,]<- c(ECCQ_stN1[i,],emosQ_N1_N2_N3_obs[[1]][[i]][18])
  ECCQ_stN2_obs[i,]<-c(ECCQ_stN2[i,],emosQ_N1_N2_N3_obs[[2]][[i]][18])
  ECCQ_stN3_obs[i,]<-c(ECCQ_stN3[i,],emosQ_N1_N2_N3_obs[[3]][[i]][18])
  
  EmosQ_stN1[i,]<-emosQ_N1_N2_N3[[1]][[i]]
  EmosQ_stN2[i,]<-emosQ_N1_N2_N3[[2]][[i]]
  EmosQ_stN3[i,]<-emosQ_N1_N2_N3[[3]][[i]]
  
  EmosQ_stN1_obs[i,]<- c(EmosQ_stN1[i,],emosQ_N1_N2_N3_obs[[1]][[i]][18])
  EmosQ_stN2_obs[i,]<-c(EmosQ_stN2[i,],emosQ_N1_N2_N3_obs[[2]][[i]][18])
  EmosQ_stN3_obs[i,]<-c(EmosQ_stN3[i,],emosQ_N1_N2_N3_obs[[3]][[i]][18])
  
}

ECCQstN1_N2_N3_obs<-list(ECCQ_stN1_obs,ECCQ_stN2_obs,ECCQ_stN3_obs)

EmosQN1_N2_N3_obs<-list(EmosQ_stN1_obs,EmosQ_stN2_obs,EmosQ_stN3_obs)


save(sample.size,n.runs,st3,ECCQstN1_N2_N3_obs,EmosQN1_N2_N3_obs,ECCstN1_N2_N3_obs,EmosN1_N2_N3_obs, file=output_file)
