setwd("/Users/sb62/Documents/Migration/Analysis_GeographicMobility_pneumoPaper/HumanMobility_RR/")
###### Script to simulate across NVT + PCV13 Fitness effect probability of staying in the home municipality, number of municipalities visited, and distance travelled . 
library(data.table)
library(ggplot2)
library(dplyr)


# load("/Users/sb62/Documents/Migration/Mobility_transmission_SA/Mobility_transmission/version2/outputs/mrcaVec.ss.opt.RData")## loads mrcavec from raw parameter overall data from 28OCT21_MunicFitSimCheck
## Functions
SeasonFunc<-function(day,r0=1.1,amp=0.03){
  
  d<-r0+amp*(sin(-pi/2+2*pi*day/365))
  return(d)}

mean=30
var=30^2
shape=mean^2/var
scale=var/mean
##INPUT TRUE MOBILITY 
load("./modelinput_data/cdr.mat.town.one.RData")
cdr.mat.town<-cdr.mat.town.one
tranmat.cdr <- cdr.mat.town

load("./ModelProjections/data/TranMatArray.234x234.RData")
TranMatArray.1<-TranMatArray.234x234
load("./modelinput_data/pairwise_geodist.town.RData") 
load("./modelinput_data/pop_municipality.2017LS.RData") 
provs <- c("Eastern Cape","Free State","Gauteng","KwaZulu-Natal","Limpopo","Mpumalanga","North West","Northern Cape","Western Cape")
tn <- rownames(cdr.mat.town)
pairwise_geodist.town <- pairwise_geodist.town[tn,tn]
load("./ModelProjections/data/densities.RData")


####Check proportion in start
mean(diag(TranMatArray.1[which(densities>500),which(densities>500),1]))
quantile(diag(TranMatArray.1[which(densities>500),which(densities>500),1]))
hist(diag(TranMatArray.1[which(densities>500),which(densities>500),1]),breaks=200)

mean(diag(TranMatArray.1[which(densities<=50),which(densities<=50),1]))
quantile(diag(TranMatArray.1[which(densities<=50),which(densities<=50),1]))
hist(diag(TranMatArray.1[which(densities<=50),which(densities<=50),1]),breaks=200)

####  ######
## Function
simFunction<-function(tranmat,specificStart=specificStart,noGens,overdis=F,R0=1,amp=0.012){
  nloc<-234
  
  if(length(specificStart)==0){
    startloc<-sample(nloc,1)}else{startloc=specificStart}
  
  genNo<-1
  whereFinal<-where<-startloc
  timeFinal<-time<-1
  nseed=1
  
  for (jj in 2:noGens){
    Reff<-SeasonFunc(time,r0=R0,amp=amp) ## With Seasonality normal
    # Reff=R0 ### Without Seasonality
    nOffspring<-rpois(nseed,Reff)
    if(overdis){nOffspring<-rnbinom(nseed,mu=Reff,size=1)}
    tmp2<-tmp<-NULL
    if(sum(nOffspring)==0)break
    for (ll in 1:nseed){
      tmp<-c(tmp,sample(nloc,nOffspring[ll],replace=T,prob=tranmat[where[ll],]))
      tmp2<-c(tmp2,time[ll]+round(rgamma(nOffspring[ll],shape=shape,scale=scale)))
    }
    # if (length(tmp2) > 20000) {print("SKIP (too long)"); next}
    if (length(tmp2) > 80000) {print("SKIP (too long)"); next}
    
    where<-tmp
    time<-tmp2
    nseed=length(where)
    whereFinal<-c(whereFinal,where)
    timeFinal<-c(timeFinal,time)
    genNo<-c(genNo,rep(jj,length(where)))
  }
  
  out.dat<-data.frame(loc=whereFinal,time=timeFinal,gen=genNo)
  dim(out.dat)
  return(out.dat)
}
R0s<-SeasonFunc(1:720,1,0.0002)*1.27
plot(R0s)


pop.dens<-c(TRUE,FALSE)
for(pd in 1:length(pop.dens)){
  rural=pop.dens[pd]
## time in days, generation of that day.
nboot=100
noGens=60
mat.boots<-list()
plot.nmuic.list <- list()
plot.dist.list <- list()
for (boot in 1:nboot) {
  times <- NULL
  repeat {
    
    

    # specificStart <- sample(234,1,prob=c(pop2019.town))
    ####choose start location based on being rural
    if(rural==TRUE){
      samp<-which(densities>=0 & densities<=50)
      # samp<-75
    
    }else{
    samp<-which(densities>=500)
    # sampe<-20
    }
    
    specificStart <- sample(samp,1)
    

    sim<-simFunction(tranmat=tranmat.cdr,specificStart=specificStart,noGens=noGens,overdis=F,R0=1,amp=0.0002) 

    print(dim(sim)[1])
    a <- sim$time[length(sim$time)]
    times <- c(times,a)
    
    print(sim$time[length(sim$time)])
    if ( sim$time[length(sim$time)]>1000  ){break} ## Overall
  }
  
  if(dim(sim)[1]>10000){next}
  
  hist(sim$time,breaks=c(100))
  simOut<-sim
  if ( max(simOut$time) > 18250) { ngen = 18250 } else( ngen=max(simOut$time))
  npairs<-which(simOut$time==ngen)[1]
  simOut.tmp2 <- simOut[1:npairs,]
  days <- sort(unique(simOut.tmp2$time))
  ts <- hist(days,c(30),plot=F)$breaks
  dmin <- ts[1:(length(ts)-1)]
  dmax <-ts[2:length(ts)]
  print(length(dmax))
  mat.tmp <- matrix(nrow=npairs,ncol=8)
  ## Weighted mean distance across all pairs 
  for(i in 1:npairs ){
    mat.tmp[i,1] <- weighted.mean(pairwise_geodist.town[simOut.tmp2[i,1],], TranMatArray.1[simOut.tmp2[i,1],,simOut.tmp2[i,3]] )
    mat.tmp[i,2] <- simOut.tmp2[i,3]
    mat.tmp[i,3] <- simOut.tmp2[i,2]
    mat.tmp[i,4] <- tn[simOut.tmp2[i,1]]
    mat.tmp[i,5] <-pop2019.town[simOut.tmp2[i,1]]
    ## N Munic per Gen
    tmp.gen <- simOut.tmp2[1:i,]
    mat.tmp[i,6] <- length(table(tmp.gen$loc))
    ## proportion home 
    start <- simOut.tmp2$loc[1]
    mat.tmp[i,7] <- length(which(simOut.tmp2[1:i,]$loc==start))/i
    mat.tmp[i,8] <- length(which(simOut.tmp2[1:i,]$loc==20))/i 
  }
  mat.tmp <- data.table(mat.tmp)
  colnames(mat.tmp) <-c("Distance","gens","days","municipality","population","NperGen","propHome","propJoBurg")
  mat.tmp$Distance<-as.numeric(mat.tmp$Distance)
  mat.tmp$gens<-as.numeric(mat.tmp$gens)
  mat.tmp$days<-as.numeric(mat.tmp$days)
  mat.tmp$population<-as.numeric(mat.tmp$population)
  mat.tmp$NperGen <- as.numeric(mat.tmp$NperGen)
  mat.tmp$propHome <- as.numeric(mat.tmp$propHome)
  mat.tmp$propJoBurg <- as.numeric(mat.tmp$propJoBurg)
  
  
  ### Mean at each generation
  maxGen=max(mat.tmp$gens,na.rm = T)
  mat <-matrix(nrow=maxGen,ncol=7)
  for ( j in 1:maxGen){
    tmp <- subset(mat.tmp, mat.tmp$gens==j )
    mat[j,1] <- mean(tmp$Distance)
    mat[j,2] <- mean(tmp$days)
    mat[j,3] <- mat.tmp$municipality[1]
    mat[j,4] <- mat.tmp$population[1]
    mat[j,5] <-  mean(tmp$NperGen)
    mat[j,6] <-  mean(tmp$propHome)
    mat[j,7] <-  mean(tmp$propJoBurg)
    
  }
  
  mat <- data.table(mat)
  colnames(mat) <-c("Distance","days","municipality","population","NperGen","propHome","propJoBurg")
  mat$gens <- c(1:maxGen)
  mat$Distance<-as.numeric(mat$Distance)
  mat$gens<-as.numeric(mat$gens)
  mat$days<-as.numeric(mat$days)
  mat$population<-as.numeric(mat$population)
  mat$NperGen<-as.numeric(mat$NperGen)
  mat$propHome<-as.numeric(mat$propHome)
  mat$propJoBurg<-as.numeric(mat$propJoBurg)
  
  mat$boot <- boot
  mat.boots[[boot]] <- mat
  print(paste0("We are at",boot,"iteration"))
}

mat.tot.rural <- rbindlist(mat.boots)




### mean at each generation
mat.tot.fin.rural <- matrix(nrow=noGens,ncol=14)
for (i in 1:noGens){
  tmp <- subset(mat.tot.rural, mat.tot.rural$gens==i)
  mat.tot.fin.rural[i,1] <- mean(tmp$NperGen)
  mn <- mean(tmp$NperGen)
  a <- sd(tmp$NperGen)
  b <- sd(tmp$NperGen)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,2] <- mn-(1.96*b)
  mat.tot.fin.rural[i,3] <- mn+(1.96*b)
  
  mat.tot.fin.rural[i,4] <- mean(tmp$Distance)
  mn <- mean(tmp$Distance)
  a <- sd(tmp$Distance)
  b <- sd(tmp$Distance)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,5] <- mn-(1.96*b)
  mat.tot.fin.rural[i,6] <- mn+(1.96*b)
  
  mat.tot.fin.rural[i,7] <- mean(tmp$propHome)
  mn <- mean(tmp$propHome)
  a <- sd(tmp$propHome)
  b <- sd(tmp$propHome)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,8] <- mn+(1.96*b)
  mat.tot.fin.rural[i,9] <- mn+(1.96*b)
  
  mat.tot.fin.rural[i,10] <- mean(tmp$propJoBurg)
  mn <- mean(tmp$propJoBurg)
  a <- sd(tmp$propJoBurg)
  b <- sd(tmp$propJoBurg)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,11] <- mn+(1.96*b)
  mat.tot.fin.rural[i,12] <- mn+(1.96*b)
  mat.tot.fin.rural[i,13] <- i
  mat.tot.fin.rural[i,14] <- mean(tmp$days)
  
  
}
mat.tot.fin.rural <-  data.table(mat.tot.fin.rural)
colnames(mat.tot.fin.rural) <- c("NperGen","NperGen_lower","NperGen_upper",
                               "Distance","Distance_lower","Distance_upper",
                               "propHome","propHome_lower","propHome_upper",
                               "propJoburg","propJoburg_lower","propJoburg_upper","gen","days")
if(rural==TRUE){ 
  mat.tot.rural.1<-mat.tot.rural
  mat.tot.fin.rural.1<-mat.tot.fin.rural
  plot.nmuic.rural <- ggplot() +
  geom_line(data=mat.tot.rural,aes(days/365, NperGen,group=as.character(boot)) ,color="grey")+
  geom_line(data=mat.tot.fin.rural.1,aes(days/365, NperGen),color="black" )+
  theme_bw()+
  # geom_ribbon(aes(ymin=lowerCI_dist,ymax=upperCI_dist))+
  scale_color_discrete(name="iter")+
  theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
  xlab("Years")+
  ggtitle("rural")+
  ylim(1,30)+
  ylab("Municipalities\n Visited (N)")
  plot.phome.rural <- ggplot() +
    geom_line(data=mat.tot.rural.1,aes(days/365, propHome,group=as.character(boot)) ,color="grey")+
    geom_line(data=mat.tot.fin.rural.1,aes(days/365, propHome),color="black" )+
    theme_bw()+
    # geom_ribbon(aes(ymin=lowerCI_dist,ymax=upperCI_dist))+
    scale_color_discrete(name="iter")+
    theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
    xlab("Years")+
    ggtitle("rural")+
    # ylim(1,30)+
    ylab("Proportion Home")
  }else{
    mat.tot.urban<-mat.tot.rural
    mat.tot.fin.urban<-mat.tot.fin.rural
    plot.nmuic.urban <- ggplot() +
      geom_line(data=mat.tot.urban,aes(days/365, NperGen,group=as.character(boot)) ,color="grey")+
      geom_line(data=mat.tot.fin.urban,aes(days/365, NperGen),color="black" )+
      theme_bw()+
      # geom_ribbon(aes(ymin=lowerCI_dist,ymax=upperCI_dist))+
      scale_color_discrete(name="iter")+
      theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
      xlab("Years")+
      ggtitle("urban")+
      ylim(1,30)+
      ylab("Municipalities\n Visited (N)")
    plot.phome.urban <- ggplot() +
      geom_line(data=mat.tot.urban,aes(days/365, propHome,group=as.character(boot)) ,color="grey")+
      geom_line(data=mat.tot.fin.urban,aes(days/365, propHome),color="black" )+
      theme_bw()+
      # geom_ribbon(aes(ymin=lowerCI_dist,ymax=upperCI_dist))+
      scale_color_discrete(name="iter")+
      theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
      xlab("Years")+
      ggtitle("urban")+
      # ylim(1,30)+
      ylab("Proportion Home")
  }
  

}
library(patchwork)
  plot.nmuic.rural+plot.nmuic.urban
  plot.phome.rural+plot.phome.urban
  
  mat.1munic<-matrix(nrow=2,ncol=3)
###percemt with 1 municipality overall
  mat.1munic[1,1]<-data.table(table(mat.tot.rural.1$NperGen==1)/nrow(mat.tot.rural.1))$N[2]
  mat.1munic[2,1]<-data.table(table(mat.tot.urban$NperGen==1)/nrow(mat.tot.urban))$N[2]
  
### percent still with 1 municipality at 2 years and 5 years
  df.rur.2<-subset(mat.tot.rural.1,mat.tot.rural.1$days>657&mat.tot.rural.1$days<803)
  mat.1munic[1,2]<-data.table(table(df.rur.2$NperGen==1)/nrow(df.rur.2))$N[2]
  df.rur.5<-subset(mat.tot.rural.1,mat.tot.rural.1$days>1642&mat.tot.rural.1$days<2008)
  mat.1munic[1,3]<-data.table(table(df.rur.5$NperGen==1)/nrow(df.rur.5))$N[2]
  
  df.urb.2<-subset(mat.tot.urban,mat.tot.urban$days>657&mat.tot.urban$days<803)
  mat.1munic[2,2]<-data.table(table(df.urb.2$NperGen==1)/nrow(df.urb.2))$N[2]
  df.urb.5<-subset(mat.tot.urban,mat.tot.urban$days>1642&mat.tot.urban$days<2008)
  mat.1munic[2,3]<-data.table(table(df.urb.5$NperGen==1)/nrow(df.urb.5))$N[2]
  
  mat.1munic<-data.table(mat.1munic)
  colnames(mat.1munic)<-c("Overall","2_years","5_years")
  mat.1munic$dens<-c("Rural","Urban")
  
  mat.1munic<-melt(mat.1munic)
  
  ggplot(mat.1munic,aes(x=variable,y=value,color=dens))+
    geom_point(size=3)+
    xlab("")+
    ylab("Proportion in 1 municipality")+
    theme_classic()+
    theme(axis.text = element_text(size=20),axis.title = element_text(size=20))
    
  
  
mat.tot.NVT <- rbindlist(mat.boots)


##########################################
######No sampling by starting location #######
###################################

## time in days, generation of that day.
nboot=100
noGens=60
mat.boots<-list()
plot.nmuic.list <- list()
plot.dist.list <- list()
for (boot in 1:nboot) {
  times <- NULL
  repeat {
    
    
    
    # specificStart <- sample(234,1,prob=c(pop2019.town))
    ####choose start location based on being rural
   
    
    specificStart <- sample(1:234,1)
    
    
    sim<-simFunction(tranmat=tranmat.cdr,specificStart=specificStart,noGens=noGens,overdis=F,R0=1,amp=0.0002) 
    
    print(dim(sim)[1])
    a <- sim$time[length(sim$time)]
    times <- c(times,a)
    
    print(sim$time[length(sim$time)])
    if ( sim$time[length(sim$time)]>1000  ){break} ## Overall
  }
  
  if(dim(sim)[1]>10000){next}
  
  hist(sim$time,breaks=c(100))
  simOut<-sim
  if ( max(simOut$time) > 18250) { ngen = 18250 } else( ngen=max(simOut$time))
  npairs<-which(simOut$time==ngen)[1]
  simOut.tmp2 <- simOut[1:npairs,]
  days <- sort(unique(simOut.tmp2$time))
  ts <- hist(days,c(30),plot=F)$breaks
  dmin <- ts[1:(length(ts)-1)]
  dmax <-ts[2:length(ts)]
  print(length(dmax))
  mat.tmp <- matrix(nrow=npairs,ncol=8)
  ## Weighted mean distance across all pairs 
  for(i in 1:npairs ){
    mat.tmp[i,1] <- weighted.mean(pairwise_geodist.town[simOut.tmp2[i,1],], TranMatArray.1[simOut.tmp2[i,1],,simOut.tmp2[i,3]] )
    mat.tmp[i,2] <- simOut.tmp2[i,3]
    mat.tmp[i,3] <- simOut.tmp2[i,2]
    mat.tmp[i,4] <- tn[simOut.tmp2[i,1]]
    mat.tmp[i,5] <-pop2019.town[simOut.tmp2[i,1]]
    ## N Munic per Gen
    tmp.gen <- simOut.tmp2[1:i,]
    mat.tmp[i,6] <- length(table(tmp.gen$loc))
    ## proportion home 
    start <- simOut.tmp2$loc[1]
    mat.tmp[i,7] <- length(which(simOut.tmp2[1:i,]$loc==start))/i
    mat.tmp[i,8] <- length(which(simOut.tmp2[1:i,]$loc==20))/i 
  }
  mat.tmp <- data.table(mat.tmp)
  colnames(mat.tmp) <-c("Distance","gens","days","municipality","population","NperGen","propHome","propJoBurg")
  mat.tmp$Distance<-as.numeric(mat.tmp$Distance)
  mat.tmp$gens<-as.numeric(mat.tmp$gens)
  mat.tmp$days<-as.numeric(mat.tmp$days)
  mat.tmp$population<-as.numeric(mat.tmp$population)
  mat.tmp$NperGen <- as.numeric(mat.tmp$NperGen)
  mat.tmp$propHome <- as.numeric(mat.tmp$propHome)
  mat.tmp$propJoBurg <- as.numeric(mat.tmp$propJoBurg)
  
  
  ### Mean at each generation
  maxGen=max(mat.tmp$gens,na.rm = T)
  mat <-matrix(nrow=maxGen,ncol=7)
  for ( j in 1:maxGen){
    tmp <- subset(mat.tmp, mat.tmp$gens==j )
    mat[j,1] <- mean(tmp$Distance)
    mat[j,2] <- mean(tmp$days)
    mat[j,3] <- mat.tmp$municipality[1]
    mat[j,4] <- mat.tmp$population[1]
    mat[j,5] <-  mean(tmp$NperGen)
    mat[j,6] <-  mean(tmp$propHome)
    mat[j,7] <-  mean(tmp$propJoBurg)
    
  }
  
  mat <- data.table(mat)
  colnames(mat) <-c("Distance","days","municipality","population","NperGen","propHome","propJoBurg")
  mat$gens <- c(1:maxGen)
  mat$Distance<-as.numeric(mat$Distance)
  mat$gens<-as.numeric(mat$gens)
  mat$days<-as.numeric(mat$days)
  mat$population<-as.numeric(mat$population)
  mat$NperGen<-as.numeric(mat$NperGen)
  mat$propHome<-as.numeric(mat$propHome)
  mat$propJoBurg<-as.numeric(mat$propJoBurg)
  
  mat$boot <- boot
  mat.boots[[boot]] <- mat
  print(paste0("We are at",boot,"iteration"))
}

mat.tot.rural <- rbindlist(mat.boots)




### mean at each generation
mat.tot.fin.rural <- matrix(nrow=noGens,ncol=14)
for (i in 1:noGens){
  tmp <- subset(mat.tot.rural, mat.tot.rural$gens==i)
  mat.tot.fin.rural[i,1] <- mean(tmp$NperGen)
  mn <- mean(tmp$NperGen)
  a <- sd(tmp$NperGen)
  b <- sd(tmp$NperGen)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,2] <- mn-(1.96*b)
  mat.tot.fin.rural[i,3] <- mn+(1.96*b)
  
  mat.tot.fin.rural[i,4] <- mean(tmp$Distance)
  mn <- mean(tmp$Distance)
  a <- sd(tmp$Distance)
  b <- sd(tmp$Distance)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,5] <- mn-(1.96*b)
  mat.tot.fin.rural[i,6] <- mn+(1.96*b)
  
  mat.tot.fin.rural[i,7] <- mean(tmp$propHome)
  mn <- mean(tmp$propHome)
  a <- sd(tmp$propHome)
  b <- sd(tmp$propHome)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,8] <- mn+(1.96*b)
  mat.tot.fin.rural[i,9] <- mn+(1.96*b)
  
  mat.tot.fin.rural[i,10] <- mean(tmp$propJoBurg)
  mn <- mean(tmp$propJoBurg)
  a <- sd(tmp$propJoBurg)
  b <- sd(tmp$propJoBurg)/sqrt(nrow(tmp))
  mat.tot.fin.rural[i,11] <- mn+(1.96*b)
  mat.tot.fin.rural[i,12] <- mn+(1.96*b)
  mat.tot.fin.rural[i,13] <- i
  mat.tot.fin.rural[i,14] <- mean(tmp$days)
  
  
}
mat.tot.fin.rural <-  data.table(mat.tot.fin.rural)
colnames(mat.tot.fin.rural) <- c("NperGen","NperGen_lower","NperGen_upper",
                                 "Distance","Distance_lower","Distance_upper",
                                 "propHome","propHome_lower","propHome_upper",
                                 "propJoburg","propJoburg_lower","propJoburg_upper","gen","days")
  mat.tot.rural.1<-mat.tot.rural
  mat.tot.fin.rural.1<-mat.tot.fin.rural
  plot.nmuic.rural <- ggplot() +
    geom_line(data=mat.tot.rural,aes(days/365, NperGen,group=as.character(boot)) ,color="grey")+
    geom_line(data=mat.tot.fin.rural.1,aes(days/365, NperGen),color="black" )+
    theme_bw()+
    # geom_ribbon(aes(ymin=lowerCI_dist,ymax=upperCI_dist))+
    scale_color_discrete(name="iter")+
    theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
    xlab("Years")+
    ggtitle("rural")+
    ylim(1,30)+
    ylab("Municipalities\n Visited (N)")
  plot.phome.rural <- ggplot() +
    geom_line(data=mat.tot.rural.1,aes(days/365, propHome,group=as.character(boot)) ,color="grey")+
    geom_line(data=mat.tot.fin.rural.1,aes(days/365, propHome),color="black" )+
    theme_bw()+
    # geom_ribbon(aes(ymin=lowerCI_dist,ymax=upperCI_dist))+
    scale_color_discrete(name="iter")+
    theme(legend.position = "none",axis.text = element_text(size=20),axis.title = element_text(size=20))+
    xlab("Years")+
    ggtitle("rural")+
    # ylim(1,30)+
    ylab("Proportion Home")
  
