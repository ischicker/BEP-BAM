rm(list=ls())
library(fields)
library(vegan)

## Minimum spanning tree ranks 
mst.rank <- function (x) {
	l.mst <- NULL
	for(f in 1:(dim(x)[2])) {
		euc.dist <- rdist(t(x[,-f]))
		l.mst <- c(l.mst,sum(spantree(euc.dist)$dist))
	}
        x.rank <- rank(l.mst,ties="random")
	return(x.rank)
}

## Multivariate ranks 
mv.rank <- function(x)
   {
     d <- dim(x)
     x.prerank <- numeric(d[2])
     for(i in 1:d[2]) {
       x.prerank[i] <- sum(apply(x<=x[,i],2,all))
     }
     x.rank <- rank(x.prerank,ties="random")
     return(x.rank)
   }

## Average ranks
avg.rank <- function(x)  {
    x.ranks <- apply(x,1,rank)
    x.preranks <- apply(x.ranks,1,mean)
    x.rank <- rank(x.preranks,ties="random")
    return(x.rank)
}

## Band depth ranks
bd.rank <- function(x)
  {
    d <- dim(x)
    x.prerank <- array(NA,dim=d)
    for(i in 1:d[1]) {
      tmp.ranks <- rank(x[i,])
      x.prerank[i,] <- (d[2] - tmp.ranks) * (tmp.ranks - 1)
    }
    x.rank <- apply(x.prerank,2,mean) + d[2] - 1
    x.rank <- rank(x.rank,ties="random")
    return(x.rank)
  } 

create.norm.data <- function(curves=19,knots=5,reps=10000,d=0,r=1)
  {
    B <- list()
    ## first row of matrix is the observation
    for(i in 1:reps){
      A <- matrix(data=rnorm(curves*knots,mean=d,sd=r),nrow=curves,ncol=knots)
      a <- rnorm(knots)
      B[[i]] <- rbind(a,A)
    }
    return(B)
  }

create.AR.data <- function(curves=19,knots=5,reps=10000,rho=1,CF.obs=function(h) exp(-h/3))
  {
    K.f <- exp(-as.matrix(dist(1:knots))/rho)  # forecasts: exponential CF with range rho
    K.o <- CF.obs(as.matrix(dist(1:knots)))  # observations: default is exponential CF with range 3
    R.f <- chol(K.f)
    R.o <- chol(K.o)

    B <- list()
    ## first row of matrix is the observation
    for(i in 1:reps){
      A <- crossprod(R.f,matrix(rnorm(curves*knots),knots,curves))
      a <- crossprod(R.o,rnorm(knots))
      B[[i]] <- t(cbind(a,A))
    }
    return(B)
  }

## Average rank histograms
avg.rhist <- function(B,reps=10000,hist_xlab="",hist_ylab="",hist_ylim=NULL)
  {
    x <- rep(0,reps) 
    for(i in 1:reps){
      B.ranks <- apply(B[[i]],2,rank)
      B.preranks <- apply(B.ranks,1,mean)
      x[i] <- rank(B.preranks,ties="random")[1]
    }
    hist(x,breaks=seq(0,20,by=1),main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray40",border="white",ylim=hist_ylim)
  }

## Multivariate rank histograms
mrh.rhist <- function(B,reps=10000,hist_xlab="",hist_ylab="",hist_ylim=NULL)
  {
    x <- rep(0,reps) 
    for(i in 1:reps){
      x[i] <- mv.rank(t(B[[i]]))[1]
    }
    hist(x,breaks=seq(0,20,by=1),main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray40",border="white",ylim=hist_ylim)
  }

## Band depth rank histograms
bd.rhist <- function(B,reps=10000,hist_xlab="",hist_ylab="",hist_ylim=NULL)
  {
    x <- rep(0,reps) 
    for(i in 1:reps){
      x[i] <- bd.rank(t(B[[i]]))[1]
    }
    hist(x,breaks=seq(0,20,by=1),main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray40",border="white",ylim=hist_ylim)
  }

## Minimum spanning tree rank histograms
mst.rhist <- function(B,reps=10000,hist_xlab="",hist_ylab="",hist_ylim=NULL)
  {
    x <- rep(0,reps) 
    for(i in 1:reps){
      x[i] <- mst.rank(t(B[[i]]))[1]
    }
    hist(x,breaks=seq(0,20,by=1),main="",xlab=hist_xlab,ylab=hist_ylab,axes=FALSE,col="gray40",border="white",ylim=hist_ylim)
  }

# ## Figures with simulated data 
# 
# pdf(file="../images/avg_rank_dim3.pdf",width=6,height=2,points=12)
# par(mfrow=c(2,3),mex=0.5,mar=c(2.7,2.7,0.1,0.1)+0.1,mgp=c(0.5,0,0))
# B <- create.norm.data(d=0,r=0.5,knots=3)
# avg.rhist(B,hist_ylab="mean = 0")
# B <- create.norm.data(d=0,r=1,knots=3)
# avg.rhist(B)
# B <- create.norm.data(d=0,r=2,knots=3)
# avg.rhist(B)
# B <- create.norm.data(d=1,r=0.5,knots=3)
# avg.rhist(B, hist_ylab="mean = 1",hist_xlab="sd = 0.5")
# B <- create.norm.data(d=1,r=1,knots=3)
# avg.rhist(B, hist_xlab="sd = 1")
# B <- create.norm.data(d=1,r=2,knots=3)
# avg.rhist(B, hist_xlab="sd = 2")
# dev.off()
# 
# pdf(file="../images/bd_rank_dim3.pdf",width=6,height=2,points=12)
# par(mfrow=c(2,3),mex=0.5,mar=c(2.7,2.7,0.1,0.1)+0.1,mgp=c(0.5,0,0))
# B <- create.norm.data(d=0,r=0.5,knots=3)
# bd.rhist(B,hist_ylab="mean = 0")
# B <- create.norm.data(d=0,r=1,knots=3)
# bd.rhist(B)
# B <- create.norm.data(d=0,r=2,knots=3)
# bd.rhist(B)
# B <- create.norm.data(d=1,r=0.5,knots=3)
# bd.rhist(B,hist_ylab="mean = 1",hist_xlab="sd = 0.5")
# B <- create.norm.data(d=1,r=1,knots=3)
# bd.rhist(B,hist_xlab="sd = 1")
# B <- create.norm.data(d=1,r=2,knots=3)
# bd.rhist(B,hist_xlab="sd = 2")
# dev.off()
# 
# png(file="../images/sd0_5.png",width=8,height=2,points=12)
# par(mfrow=c(2,4),mex=0.5,mar=c(2.7,2.7,0.1,0.1)+0.1,mgp=c(0.5,0,0))
B <- create.norm.data(d=0,knots=5,r=0.5)
# mrh.rhist(B,hist_ylab="d = 5")
# avg.rhist(B)
# bd.rhist(B)
# mst.rhist(B)
# B <- create.norm.data(d=0,knots=15,r=0.5)
# mrh.rhist(B,hist_ylab="d = 15", hist_xlab="Multivariate Rank")
# avg.rhist(B,hist_xlab="Average Rank")
# bd.rhist(B,hist_xlab="Band Depth Rank")
# mst.rhist(B,hist_xlab="Minimum Spanning Tree Rank")
# dev.off()
# 
# pdf(file="../images/sd2.pdf",width=8,height=2,points=12)
# par(mfrow=c(2,4),mex=0.5,mar=c(2.7,2.7,0.1,0.1)+0.1,mgp=c(0.5,0,0))
# B <- create.norm.data(d=0,knots=5,r=2)
# mrh.rhist(B,hist_ylab="d = 5")
# avg.rhist(B)
# bd.rhist(B)
# mst.rhist(B)
# B <- create.norm.data(d=0,knots=15,r=2)
# mrh.rhist(B,hist_ylab="d = 15", hist_xlab="Multivariate Rank")
# avg.rhist(B,hist_xlab="Average Rank")
# bd.rhist(B,hist_xlab="Band Depth Rank")
# mst.rhist(B,hist_xlab="Minimum Spanning Tree Rank")
# dev.off()
# 
# pdf(file="../images/AR1_dim5.pdf",width=8,height=2.5,points=12)
# par(mfrow=c(2,4),mex=0.5,mar=c(2.7,2.7,0.1,0.1)+0.1,mgp=c(0.5,0,0))
# B <- create.AR.data(rho=2)
# y.max <- 700
# mrh.rhist(B,hist_ylab=expression(paste(tau, " = 2")),hist_ylim=c(0,y.max))
# avg.rhist(B,hist_ylim=c(0,y.max))
# bd.rhist(B,hist_ylim=c(0,y.max))
# mst.rhist(B,hist_ylim=c(0,y.max))
# B <- create.AR.data(rho=4)
# mrh.rhist(B,hist_ylim=c(0,y.max),hist_ylab=expression(paste(tau," = 4")),hist_xlab="Multivariate Rank")
# avg.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Average Rank")
# bd.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Band Depth Rank")
# mst.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Minimum Spanning Tree Rank")
# dev.off()
# 
# pdf(file="../images/AR1vsOthers.pdf",width=8,height=3.7,points=12)
# par(mfrow=c(3,4),mex=0.5,mar=c(2.7,2.7,0.1,0.1)+0.1,mgp=c(0.5,0,0))
# y.max <- 850
# B <- create.AR.data(knots=15,rho=3,CF.obs=function(h) exp(-h/4.5)*(0.75+0.25*cos(h*pi/2)))
# mrh.rhist(B,hist_ylim=c(0,y.max),hist_ylab="corr. model a)",hist_xlab="Multivariate Rank")
# avg.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Average Rank")
# bd.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Band Depth Rank")
# mst.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Minimum Spanning Tree Rank")
# B <- create.AR.data(knots=15,rho=3,CF.obs=function(h) 1/(1+h/2.5))
# mrh.rhist(B,hist_ylim=c(0,y.max),hist_ylab="corr. model b)",hist_xlab="Multivariate Rank")
# avg.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Average Rank")
# bd.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Band Depth Rank")
# mst.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Minimum Spanning Tree Rank")
# B <- create.AR.data(knots=15,rho=3,CF.obs=function(h) pmax(1-h/5,0))
# mrh.rhist(B,hist_ylim=c(0,y.max),hist_ylab="corr. model c)",hist_xlab="Multivariate Rank")
# avg.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Average Rank")
# bd.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Band Depth Rank")
# mst.rhist(B,hist_ylim=c(0,y.max),hist_xlab="Minimum Spanning Tree Rank")
# dev.off()

## Tables 1 and 2

## Average ranking
avg.stats <- function(curves=19,knots=5,reps=30000,rho=1,hist_xlab="",hist_ylab="")
{
	x <- rep(0,reps)
	z <- rep(0,reps)
	x.pre <- rep(0,reps)
	z.pre <- rep(0,reps)
	
	K.f <- exp(-as.matrix(dist(1:knots))/rho)  # forecasts: exponential CF with range rho
	K.o <- exp(-as.matrix(dist(1:knots))/3)  # observations: exponential CF range 1
	R.f <- chol(K.f)
	R.o <- chol(K.o)
	
## first row of matrix is the observation
	for(i in 1:reps){
		A <- crossprod(R.f,matrix(rnorm(curves*knots),knots,curves))
		a <- crossprod(R.o,rnorm(knots))
		B <- cbind(a,A)
		
		B.ranks <- apply(B,1,rank)
		B.preranks <- apply(B.ranks,1,mean)
		x.pre[i] <- B.preranks[1]
		z.pre[i] <- B.preranks[2]
		
		tmp <- rank(B.preranks,ties="random")
		x[i] <- tmp[1]
		z[i] <- tmp[2]
	}
	print("Average ranks")
	print(paste("Observation: ", round(mean(x), 1), round(var(x)) ))
	print(paste("Forecast: ", round(mean(z), 1), round(var(z)) ))
	
}

## Band depth ranking
bd.stats <- function(curves=19,knots=5,reps=30000,rho=1,hist_xlab="",hist_ylab="")
{
	x <- rep(0,reps)
	z <- rep(0, reps)
	x.pre <- rep(0,reps)
	z.pre <- rep(0,reps)
	
	K.f <- exp(-as.matrix(dist(1:knots))/rho)  # forecasts: exponential CF with range rho
	K.o <- exp(-as.matrix(dist(1:knots))/3)    # observations: exponential CF range 1
	R.f <- chol(K.f)
	R.o <- chol(K.o)
	
## first row of matrix is the observation
	for(i in 1:reps){
		A <- crossprod(R.f,matrix(rnorm(curves*knots),knots,curves))
		a <- crossprod(R.o,rnorm(knots))
		B <- cbind(a,A)
		
		tmp <- bd.rank(B)
		x[i] <- tmp$x.rank[1]
		z[i] <- tmp$x.rank[2]
		x.pre[i] <- tmp$x.prerank[1]
		z.pre[i] <- tmp$x.prerank[2]
		
	}
	
	print("Band depth ranks")
	print(paste("Observation: ", round(mean(x), 1), round(var(x)) ))
	print(paste("Forecast: ", round(mean(z), 1), round(var(z)) ))
}

# print("rho=2, m=20, d=5")
# avg.stats(rho=2)
# bd.stats(rho=2)
# print("rho=2, m=100, d=5")
# avg.stats(rho=2, curves=99)
# bd.stats(rho=2, curves=99)
# print("rho=2, m=200, d=5")
# avg.stats(rho=2, curves=199)
# bd.stats(rho=2, curves=199)
# print("rho=2, m=500, d=5")
# avg.stats(rho=2, curves=499)
# bd.stats(rho=2, curves=499)
# print("rho=2, m=20, d=100")
# avg.stats(rho=2, knots=100)
# bd.stats(rho=2, knots=100)
# print("rho=2, m=100, d=100")
# avg.stats(rho=2, curves=99, knots=100)
# bd.stats(rho=2, curves=99, knots=100)
# print("rho=2, m=200, d=100")
# avg.stats(rho=2, curves=199, knots=100)
# bd.stats(rho=2, curves=199, knots=100)
# print("rho=2, m=500, d=100")
# avg.stats(rho=2, curves=499, knots=100)
# bd.stats(rho=2, curves=499, knots=100)
# print("rho=2, m=20, d=200")
# avg.stats(rho=2, knots=200)
# bd.stats(rho=2, knots=200)
# print("rho=2, m=100, d=200")
# avg.stats(rho=2, curves=99, knots=200)
# bd.stats(rho=2, curves=99, knots=200)
# print("rho=2, m=200, d=200")
# avg.stats(rho=2, curves=199, knots=200)
# bd.stats(rho=2, curves=199, knots=200)
# print("rho=2, m=500, d=200")
# avg.stats(rho=2, curves=499, knots=200)
# bd.stats(rho=2, curves=499, knots=200)
# print("rho=2, m=20, d=500")
# avg.stats(rho=2, knots=500)
# bd.stats(rho=2, knots=500)
# print("rho=2, m=100, d=500")
# avg.stats(rho=2, curves=99, knots=500)
# bd.stats(rho=2, curves=99, knots=500)
# print("rho=2, m=200, d=500")
# avg.stats(rho=2, curves=199, knots=500)
# bd.stats(rho=2, curves=199, knots=500)
# print("rho=2, m=500, d=500")
# avg.stats(rho=2, curves=499, knots=500)
# bd.stats(rho=2, curves=499, knots=500)

## Postprocessing of ECMWF weather forecasts 
load("temp.obs.Rdata") 

stat.id <- "10382"  

first.day <- 20101001
last.day  <- 20121231

lead.time <- seq(6,72,6)   

td <- 50    # training days
es <- 50    # ensemble size

## The ECMWF ensemble forecasts are freely availble from the TIGGE repository, 
## see https://protect-us.mimecast.com/s/ANrmBzULqR9rc9
## We obtain a station-specific forecast by a bilinear interpolation of the 
## forecasts at the four nearest grid points. 
## The SYNOP station Berlin Tegel is located at latitude 52.5644 and  
## longitude 13.3088. 
load(paste("station",stat.id,".Rdata",sep=""))

forecast.index <- which( days >= first.day & days <= last.day )
station.index <- which( temp.obs$station == stat.id )

l <- length(forecast.index)
d <- length(lead.time)

training.table <- vector(l, mode="list")
names(training.table) <- days[forecast.index]

for (i in 1:l)  {
	training.table[[i]] <- tail( days[days<days[forecast.index][i]], 61 )
}

ens.mean <- apply(ens[,,1+lead.time%/%3]-273.15, c(1,3), mean, na.rm=TRUE)

# create an observation matrix of the same size as 'ens.mean' and 'ens.var', taking the lead
#  time into account, i.e. the i-th row and j-th column corresponds to the forecast initialized
#  at 0000 UTC on days[i] with lead.time[j]

obs <- matrix(NA,dim(ens)[1],d)
for (j in 1:d)  {
	hh.index <- which( temp.obs$date %% 100 == lead.time[j] %% 24 )
	obs.hh.stat <- temp.obs$obs[hh.index,station.index]
	obs.date <- temp.obs$date[hh.index] %/% 100
	shift <- lead.time[j] %/% 24
	obs[,j] <- obs.hh.stat[ which(obs.date >= days[1] & obs.date <= tail(days,1)) + shift ]
}

####   Compute bias correction, create error dressing ensemble 
biasCoef <- array(dim=c(l,d,2))
forcMean <- matrix(NA,l,d)

ED.ens.ind <- ED.ens.ecc <- ED.ens.vec.mnorm <- array(dim=c(l,es,d))

errors <- matrix(NA,d,td)

set.seed(17)

for (i in 1:l)
{
    verif.day <- days[forecast.index][i]
    day.index <- which(days==verif.day)
    
    cat(paste("\n\nModel fitting for date", verif.day, "\n"))
    sample.ind <- sample(1:td, es, replace=TRUE)
    
    tmp <- rep(0,d)
    for (j in 1:d)
	{
        shift <- max(lead.time) %/% 24
        training.index <- which( days %in% training.table[[i]][(-(td-1):0)+length(training.table[[i]])-shift] )
        
        est <- lm ( obs[training.index,j] ~ ens.mean[training.index,j] )
        
        biasCoef[i,j,] <- est$coef
        forcMean[i,j] <- biasCoef[i,j,1] + biasCoef[i,j,2]*ens.mean[day.index,j]
        
        errors[j,] <- est$residuals
        
        X <- cbind(1,ens.mean[training.index,j])
        x0 <- c(1,ens.mean[day.index,j])
        corr.fct <- (td/(td-2))*(1 + crossprod(x0, solve(crossprod(X,X),x0)) )
        tmp[j] <- sqrt(corr.fct)
	## indep. error dressing ensemble
        ED.ens.ind[i,,j] <- forcMean[i,j] + sqrt(corr.fct)*sample(errors[j,], es, replace=TRUE)
	## ECC error dressing ensemble
        ED.ens.ecc[i,order(ens[forecast.index[i],,(1+lead.time%/%3)[j]]),j] <- sort(ED.ens.ind[i,,j])
	}
    
    for(k in 1:50)
	{
	## multivariate normal model
        ED.ens.vec.mnorm[i,k,] <- forcMean[i,] + tmp * (rnorm(d) %*% chol(cov(t(errors))))
	}
}

#### Univariate Histograms
lead.range <- 1+range((lead.time-1) %/% 24)
vr.ED <- 1 + apply( sweep(ED.ens.ind,c(1,3),obs[forecast.index,],"<="), c(1,3), sum)

pdf(file=paste("../images/ED-univ-",stat.id,"-day",lead.range[1],"-day",lead.range[2], ".pdf",sep=""),width=8,height=3.5,points=12)
par(mfrow=c(3,4),mex=0.5,mar=c(2.7,0.1,0.1,0.1)+0.1,mgp=c(0.5,0,0))
for (j in 1:d) {
	hist(vr.ED[,j],breaks=seq(0.5,51.5,1),main="",xlab=paste("Lead Time = ",lead.time[j],"h",sep=""),ylab="",axes=FALSE,col="gray40",border="white",ylim=c(0,35))
}
dev.off()


####  Multivariate Histograms
mv.rank.ED.ind <- mbd.rank.ED.ind <- avg.rank.ED.ind <- mst.rank.ED.ind <- numeric(l)
mv.rank.ED.ecc <- mbd.rank.ED.ecc <- avg.rank.ED.ecc <- mst.rank.ED.ecc <- numeric(l)
mv.rank.ED.vec.mn <- mbd.rank.ED.vec.mn <- avg.rank.ED.vec.mn <- mst.rank.ED.vec.mn <- numeric(l)

for (i in 1:l)  {
	data <- cbind(obs[forecast.index[i],],t(ED.ens.ind[i,,]))
	mv.rank.ED.ind[i] <- mv.rank(data)[1]
	mbd.rank.ED.ind[i] <- bd.rank(data)[1]
	avg.rank.ED.ind[i] <- avg.rank(data)[1]
	mst.rank.ED.ind[i] <- mst.rank(data)[1]
}

for (i in 1:l)  {
	data <- cbind(obs[forecast.index[i],],t(ED.ens.ecc[i,,]))
	mv.rank.ED.ecc[i] <- mv.rank(data)[1]
	mbd.rank.ED.ecc[i] <- bd.rank(data)[1]
	avg.rank.ED.ecc[i] <- avg.rank(data)[1]
	mst.rank.ED.ecc[i] <- mst.rank(data)[1]
}

for (i in 1:l)  {
	data <- cbind(obs[forecast.index[i],],t(ED.ens.vec.mnorm[i,,]))
	mv.rank.ED.vec.mn[i] <- mv.rank(data)[1]
	mbd.rank.ED.vec.mn[i] <- bd.rank(data)[1]
	avg.rank.ED.vec.mn[i] <- avg.rank(data)[1]
	mst.rank.ED.vec.mn[i] <- mst.rank(data)[1]
}

pdf(file="../images/ED-BerlinTegel-multi2.pdf", width=8, height=4,points=12)
par(mfrow=c(3,4),mex=0.5,mar=c(2.5,2.2,0,0)+0.1,mgp=c(0.5,0,0))
K <- 125
hist(mv.rank.ED.ind,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="Independent",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(avg.rank.ED.ind,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mbd.rank.ED.ind,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mst.rank.ED.ind,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mv.rank.ED.ecc,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="ECC",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(avg.rank.ED.ecc,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mbd.rank.ED.ecc,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mst.rank.ED.ecc,breaks=seq(0.5,51.5,1),main="",xlab="",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mv.rank.ED.vec.mn,breaks=seq(0.5,51.5,1),main="",xlab="Multivariate Rank",ylab="Multivariate Normal",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(avg.rank.ED.vec.mn,breaks=seq(0.5,51.5,1),main="",xlab="Average Rank",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mbd.rank.ED.vec.mn,breaks=seq(0.5,51.5,1),main="",xlab="Band Depth Rank",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
hist(mst.rank.ED.vec.mn,breaks=seq(0.5,51.5,1),main="",xlab="Minimum Spanning Tree Rank",ylab="",axes=FALSE,col="gray40",border="white", ylim=c(0,K))
dev.off()







