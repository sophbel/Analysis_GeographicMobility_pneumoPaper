####Fit model with simulation
#########Set up simulation
library(dplyr)
library(data.table)
library(ggplot2)
library(abind)
library(doParallel)
library(ucminf)
library(doMC)
library(Rcpp)
library(RcppEigen)
library(Rfast)
library(coda)
library(fmcmc)
setwd("/Users/sb62/Documents/Migration/SA_Migration_110422/RData_outputs/")
load('cdr.mat.one.RData')
cdr.mat<-cdr.mat.one
load('cdr.mat.town.one.RData')
cdr.mat.town<-cdr.mat.town.one
load('pairwise_geodist.RData')
load("dat.tmp.allser.RData")
load('pop_2019.RData')
load("pop2019_municipality.2017LS.RData")
sourceCpp("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/MatrixMultiplication.cpp")

# load("/Users/sb62/Documents/Migration/Mobility_transmission_SA/Mobility_transmission/MCMC_ReviseModel/outputs/Munic_9X9/gen30/ans_munic9_2000.100gens.30genTime.RData")
# load("/Users/sb62/Documents/Migration/Mobility_transmission_SA/Mobility_transmission/MCMC_ReviseModel/outputs/Munic_9X9/GenerationTimeSensitivity/gen30/ans_munic9_3000.100gen.genTime30.RData")
# ans<-mcmc(ans.munic,start=1000, end=2000)
# fit.par<-summary(ans)$statistics[,1]

#####index for province and municipality match
nloc.munic<-length(pop2019.town)
splitnames <- matrix(nrow=nloc.munic,ncol=2)
for (us in 1:nloc.munic) {
  for (prov in 1:2) { 
    splitnames[us,prov] <- strsplit(names(pop2019.town),"_")[[us]][prov]
  }
}
Eigencpp=TRUE
a<-match(colnames(cdr.mat),names(pop_2019))
pop_2019<-pop_2019[a]

a<-match(colnames(cdr.mat.town),names(pop2019.town))
pop2019.town<-pop2019.town[a]


## Simulatation parameters
nsims<-1  ## Number of sepaarte simulations to run
burnin<-20
nGen2<-nGen<-1:30
pGen2<-pGen<-rep(1/30,30)
ntimes<-50000 ## No. of pairs per simulation pre observation

### True sampling from true data set.
dat.inMaster<-do.call("rbind",dat.tmp.allser)
probByloc<-hist(c(dat.inMaster[,3],dat.inMaster[,4]),breaks=1:10-0.5,plot=F)$counts
probSample<-probByloc
# probSample<-probSample/max(probSample)
probSample<-probSample/sum(probSample)

##### Add parameter to Mobility data 
parHome<--2
min.range=-0.04
max.range=0.6
nlocs=nloc=9
npop<-pop_2019
pHome<-npop/sum(npop)
#########
tmpbase<-cdr.mat
tmppar1 <- exp(parHome)/(1+exp(parHome))
tmppar <- min.range+tmppar1*(max.range-min.range)
tmpdiag<-diag(tmpbase)-tmppar
tmpdiag[which(tmpdiag>0.99999)]<-0.99999
diag(tmpbase)<-0
tmpbase<-sweep(tmpbase,1,rowSums(tmpbase),"/")
tmpbase<-sweep(tmpbase,1,(1-tmpdiag)/(1-diag(tmpbase)),"*")
diag(tmpbase)<-tmpdiag
tmp.sick<-tmp<-tmpbase

move3<-tcrossprod(tmp.sick,tmp.sick)
move4<-sweep(move3,2,pHome,"*")
TranMat<-sweep(move4,1,rowSums(move4),"/")


############
#SIMULATION FUNCTION
SimulationFunction<-function(){
  
  start=sample(nloc,ntimes,prob=pHome,replace=T)
  for (jj in 1:burnin){
    for (kk in 1:ntimes){
      start[kk]=sample(nloc,1,prob=TranMat[start[kk],])
    }
  }
  whereEnd2<-whereEnd<-start
  nogens<-sample(nGen,ntimes,replace=T,prob=pGen)
  nogens2<-sample(nGen2,ntimes,replace=T,prob=pGen2)
  for (kk in 1:ntimes){
    for (jj in 1:nogens[kk]){
      whereEnd[kk]=sample(nloc,1,prob=TranMat[whereEnd[kk],])
    }
    for (jj in 1:nogens2[kk]){
      whereEnd2[kk]=sample(nloc,1,prob=TranMat[whereEnd2[kk],])
    }
  }
  
  output<-list(start, whereEnd,whereEnd2,nogens,nogens2)
  return(output)
}


sim<-SimulationFunction()

start<-sim[[1]]
whereEnd<-sim[[2]]
whereEnd2<-sim[[3]]
nogens<-sim[[4]]
nogens2<-sim[[5]]


##Prepare data
dat.in.all<-as.data.frame(cbind(id1=1:length(start),id2=1:length(start),
                                start=start,loc1=whereEnd2,loc2=whereEnd,
                                nogens=nogens,nogens2=nogens2,meanGen=nogens+nogens2))

a<-as.matrix(dat.in.all[,4:5])
dat.in.all$dist<-pairwise_geodist[a]
a<-as.matrix(dat.in.all[,c(3,5)])
dat.in.all$distStart1<-pairwise_geodist[a]
a<-as.matrix(dat.in.all[,c(3,4)])
dat.in.all$distStart2<-pairwise_geodist[a]

####SUB-Sample those in the same place
samp1<-rbinom(nrow(dat.in.all)*3,1,prob=probSample[rep(dat.in.all[,"loc1"],3)])
samp2<-rbinom(nrow(dat.in.all)*3,1,prob=probSample[rep(dat.in.all[,"loc2"],3)])
# samp1<-rbinom(nrow(dat.in.all),1,prob=probSample[dat.in.all[,"loc1"]])
# samp2<-rbinom(nrow(dat.in.all),1,prob=probSample[dat.in.all[,"loc2"]])
dat.in.tmp<-dat.in.all[which(samp1==1&samp2==1),]

print(dim(dat.in.tmp))
locs1<-table(dat.in.tmp$loc1)
locs2<-table(dat.in.tmp$loc2)
print(locs1)
print(locs2)


###################
##############REFIT MODEL USING SIMULATION##########################################
calcAllProbs=FALSE
dat.in2<-dat.in.tmp

ncore=1
# min.range<-(-0.04)
# max.range<-0.6

extTranMatDat.tmp<-list()
extTranMatDat.tmp$popbyCell<-pop_2019
extTranMatDat.tmp$pars<-list()
extTranMatDat.tmp$pars$homeSus<-999

genTime<-30

varGen<-(genTime)^2
shape=genTime^2/varGen
scale=varGen/genTime

max.no.gens<-30 ### 8.2 years for minimum number of generations times (10 days)
maxGen<-max.no.gens

#### Plot generation time distribution. 
hist(rgamma(5000, shape=shape, scale=scale),breaks=500)
# 
probByGen.tmp<-array(0,c(nrow(dat.in2),maxGen,2))
for(k in 1:nrow(dat.in2)){
  probByGen.tmp[k,dat.in2[k,6],1]<-1
  probByGen.tmp[k,dat.in2[k,7],2]<-1
}


source("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/LikFunc.mcmc.SIMS.Munic.R")

# nInfecLoc<-colSums(rbind(table(dat.in.all$loc1),table(dat.in.all$loc2)))
nInfecLoc<-table(dat.in.all$start)

dat.in2<-dat.in.tmp
npairs=nrow(dat.in2)
##############################################################################################################
##RUN LIKELIHOOD FUNCTION
#####RUN LLMETHOD#######################
# singPar=FALSE#####If I fit only the diagonal or if there are the additional infection parameters included.
# if(singPar==TRUE){
# source("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/LikFunc.mcmc.SIMS.Munic.R")
# 
# pars<-seq(-2.8,-1.3,0.1)
# par.ll<-rep(NA,length(pars))
# for (i in 1:length(pars)){
# par.ll[i]<-likFunc.sim(pars[i])
# print(i)
# }
# 
# plot(pars,par.ll)
# }else
# 
# #######RUN MCMC#######################
# {
#   source("./MCMC_model/Simulations/LikFunc.mcmc.SIMS.Munic.R")
# 
# iters=5000
# startPar<-c(-2.5,rep(0,8))
# load("./MCMC_model/Simulations/output/sim.ans.RData")
# ans<-mcmc(sim.ans,start=2000,end=4000)
# rm(sim.ans)
# startPar<-summary(ans)$statistics[,1]
# start.time<-Sys.time()
# # sim.ans<-MCMC(likFunc.sim,initial = startPar,nsteps  = iters,kernel  = kernel_normal(scale = .06),progress = interactive())
# end.time<-Sys.time()
# print(end.time-start.time)
# }


#################################################################################################################################################

#####################PLOT OUTPUT 
singPar=FALSE
################Single Parameter
if(singPar==TRUE){
parHome<-pars[which(par.ll==max(par.ll))]
Timehome=exp(parHome)/(1+exp(parHome))
extTranMatDat.tmp$pars$homeSus<-parHome}
# else{###beginning of else for single par vs multiple par
################Multiple Parameter
library(fmcmc)

  # load("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/output/sim.ans.0.05.10000.RData")
  # load("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/output/sim.ans.0.06.10000.RData")
  load("./MCMC_model/Simulations/output/sim.ans.0.06.20000minmaxBackin.RData")
  # load("/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/output/sim.ans.0.06.20000.RData")
  
  plot(sim.ans[,1])
  ans<-mcmc(sim.ans,start=15000,end=20000)
  plot(ans[,1])
  summary(ans)$statistics[,1]
  1-rejectionRate(ans)
  posteriors<-as.matrix(ans)
  nsim<-2000
  pars<-matrix(nrow=nsim,ncol=ncol(ans))
  for(i in 1:nsim){
    simno<-sample(nrow(posteriors),1)
    pars[i,]<-posteriors[simno,]
  }
true.modfit<-sub.modfit<-matrix(ncol=nsim,nrow=maxGen)
for(post in 1:nsim){
par<-pars[post,1]
# par<-summary(ans)$statistics[,1]
parHome<-par[1]
extTranMatDat.tmp$pars$homeSus<-parHome
Timehome=exp(parHome)/(1+exp(parHome))
nInfecLoc<-rep(1,9)
nInfecLoc[1:8]<-exp(par[2:9])
nInfecLoc<-nInfecLoc/sum(nInfecLoc)
###No. detected
avNoDetectByloc2<-colSums(rbind(hist(dat.in2$loc1,plot=FALSE,breaks=seq(0.5,9.5,1))$counts,hist(dat.in2$loc2,plot=FALSE,breaks=seq(0.5,9.5,1))$counts))
# avNoDetectByloc2<-as.numeric(table(c(dat.in2$loc1,dat.in2$loc2)))
no.detect<-avNoDetectByloc2
probByloc<-no.detect/nInfecLoc
probByloc<-probByloc/sum(probByloc)

#####Truth
nInfecLoc.true<-colSums(rbind(hist(dat.in.all$loc1,plot=FALSE,breaks=seq(0.5,9.5,1))$counts,hist(dat.in.all$loc2,plot=FALSE,breaks=seq(0.5,9.5,1))$counts))
nInfecLoc.true<-nInfecLoc.true/sum(nInfecLoc.true)
incidence.true<-nInfecLoc.true/extTranMatDat.tmp$popbyCell
#####Check how it looks against population size
incidence<-nInfecLoc/extTranMatDat.tmp$popbyCell
tmp1<-cbind(exp(posteriors[,2:9]),1)
tmp<-sweep(tmp1,1,rowSums(tmp1),"/")
quantiles.probInfec<-apply(tmp,2,quantile,probs=c(0.025,0.5,0.975))


####Plot infection probability against the truth
# pdf(file="/Users/sb62/Documents/Migration/SA_Migration_110422/MCMC_model/Simulations/plots/infecProbLoc.pdf",height=2,width=2)
# # quartz(width=2,height=2)
# par(mar=c(3,3,1,1))
# plot(extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[2,],pch=20,axes=F,ylim=c(0,0.5))
# points(extTranMatDat.tmp$popbyCell/1000000,nInfecLoc.true,col="red",pch=20,axes=F,ylim=c(0,0.5))
# arrows(extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[1,],extTranMatDat.tmp$popbyCell/1000000,quantiles.probInfec[3,],length=0)
# axis(1,cex=0.75,tck=-0.04,padj=-1.5,cex.axis=0.75)
# axis(2,cex=0.75,tck=-0.04,hadj=0.5,cex.axis=0.75,las=1)
# mtext(2,text="Proportion of infections",cex=0.75,line=1.7)
# mtext(1,text="Population size (x10^6)",cex=0.75,line=1.2)
# dev.off()

# }####end of ifelse for single par vs. multiple par
################################

nlocs=nloc=9
npop<-pop_2019
pHome<-npop/sum(npop)
# min.range=-0.04
# max.range=0.6
# tmp.pHome.MUNIC<-pop2019.town
# tmp.pHome.MUNIC<-tmp.pHome.MUNIC/sum(tmp.pHome.MUNIC)
tmp.pHome<-extTranMatDat.tmp$popbyCell
tmp.pHome<-tmp.pHome/sum(tmp.pHome)

#### Create mobility matrix for susceptible individuals
tmpbase<-cdr.mat
tmppar1 <- Timehome
tmppar <- min.range+tmppar1*(max.range-min.range)
tmpdiag<-diag(tmpbase)-tmppar
tmpdiag[which(tmpdiag>0.99999)]<-0.99999
diag(tmpbase)<-0
tmpbase<-sweep(tmpbase,1,rowSums(tmpbase),"/")
tmpbase<-sweep(tmpbase,1,(1-tmpdiag)/(1-diag(tmpbase)),"*")
diag(tmpbase)<-tmpdiag
tmp.sick<-tmp<-tmpbase

# probInfec.tmp<-rep(1,nloc)
# tmp<-sweep(tmp,2,probInfec.tmp,"*")
move3<-tcrossprod(tmp.sick,tmp.sick)
move4<-sweep(move3,2,tmp.pHome,"*")
TranMat.tmp2<-sweep(move4,1,rowSums(move4),"/")

maxGen=30
probByGenA<-probByGenB<-rep(1,maxGen)
gensB<-gensA<-(1:maxGen)
maxgen.tmp<-max(c(gensA,gensB))
TranMatArray<-array(NA,c(nlocs,nlocs,maxgen.tmp))
TranMatArray[,,1]<-TranMat.tmp2
if(Eigencpp){
  for (j in 2:maxgen.tmp){TranMatArray[,,j]<-eigenMapMatMult(TranMatArray[,,j-1], TranMat.tmp2)
  }
}else{
  for (j in 2:maxgen.tmp){TranMatArray[,,j]<-TranMatArray[,,j-1]%*% TranMat.tmp2}
}

mrcaVec<-pHome
mrcaVecTRUE<-hist(dat.in.all$start,breaks=(0:nloc)+0.5,plot=F)$counts
mrcaVecTRUE<-mrcaVecTRUE/sum(mrcaVecTRUE)
# mrcaVec<-mrcaVecTRUE	


TranMatArrayA.PreSamp<-TranMatArrayA<-TranMatArray[,,gensA]
if(singPar==TRUE){
TranMatArrayA<-sweep(TranMatArrayA,2,probSample,"*")}else{
TranMatArrayA<-sweep(TranMatArrayA,2,probByloc,"*")}####MULTIPLE PARAMETERS


TranMatArrayB<-TranMatArray[,,gensB]
if(singPar==TRUE){
TranMatArrayB<-sweep(TranMatArrayB,2,probSample,"*")}else{
TranMatArrayB<-sweep(TranMatArrayB,2,probByloc,"*")}#####MULTIPLE PARAMETERS
TranMatArrayB.2<-sweep(TranMatArrayB,1,mrcaVec,"*")

TranMatArray.PreSamp<-sweep(TranMatArrayA.PreSamp,1,mrcaVec,"*")
TranMatArrayA2<-sweep(TranMatArrayA,1,mrcaVec,"*")
TranMatArrayA2.mrcaTRUE<-sweep(TranMatArrayA,1,mrcaVecTRUE,"*")
TranMatArrayA3<-matrix(TranMatArrayA2,nlocs,nlocs*length(gensA))
TranMatArrayA3.mrcaTRUE<-matrix(TranMatArrayA2.mrcaTRUE,nlocs,nlocs*length(gensA))
TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nlocs,byrow=T)
probAllPrs<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3)
probAllPrs.mrcaTRUE<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3.mrcaTRUE)

nogens<-rep(1:maxGen,each=nlocs)
reflocs<-rep(1:nlocs,maxGen)
nogens.mat<-outer(nogens,nogens,"+")
a<-as.matrix(expand.grid(reflocs,reflocs))
refDist<-pairwise_geodist[a]
refDist.mat<-matrix(refDist,maxGen*nlocs,maxGen*nlocs)

possGens<-expand.grid(1:maxGen,1:maxGen)
possGenstot<-possGens[,1]+possGens[,2]

distgenfromMRCA2<-distgenfromMRCA<-distgenfromMRCA.mrcaTRUE<-rep(NaN,maxGen)
for (i in 1:maxGen){
  distgenfromMRCA[i]<-weighted.mean(pairwise_geodist,TranMatArrayA2[,,i])
  distgenfromMRCA.mrcaTRUE[i]<-weighted.mean(pairwise_geodist,TranMatArrayA2.mrcaTRUE[,,i])
  distgenfromMRCA2[i]<-weighted.mean(pairwise_geodist,TranMatArray.PreSamp[,,i])
}

weighted.distgen<-rep(NaN,maxGen)
for (i in 1:maxGen){
  a<-which(nogens.mat==i)
  weighted.distgen[i]<-weighted.mean(refDist.mat[a],probAllPrs[a])
}

weighted.distgen2.mrcaTRUE<-weighted.distgen2<-rep(NaN,maxGen)
for(i in 2:maxGen){
  a<-which(possGenstot==i)
  tmpWt2<-tmpWt<-rep(NaN,length(a))
  
  for (j in 1:length(a)){
    indA<-possGens[a[j],1]
    indB<-possGens[a[j],2]
    
    probDestA<-TranMatArray[,,indA]
    probDestB<-TranMatArray[,,indB]
    

    
    probDestA<-sweep(probDestA,2,probSample,"*") ##ProbDetectedA
    probDestB<-sweep(probDestB,2,probSample,"*")
    
    # probDestA<-sweep(probDestA,2,probSample,"*") ##MULTIPLE PARS
    # probDestB<-sweep(probDestB,2,probSample,"*")
    
    
    dist.tmp<-rep(NaN,9)
    wtDist2<-wtDist<-matrix(NaN,0,2)
    for(k in 1:9){
      b<-outer(probDestA[k,],probDestB[k,],"*")*mrcaVec[k]
      # b<-b/sum(b)
      # dist.tmp[k]<-weighted.mean(pairwise_geodist,b)
      wtDist<-rbind(wtDist,cbind(as.vector(pairwise_geodist),as.vector(b)))
      b<-outer(probDestA[k,],probDestB[k,],"*")*mrcaVecTRUE[k]
      wtDist2<-rbind(wtDist2,cbind(as.vector(pairwise_geodist),as.vector(b)))
    }
    # tmpWt[j]<-mean(dist.tmp)
    tmpWt[j]<-weighted.mean(wtDist[,1],wtDist[,2])
    tmpWt2[j]<-weighted.mean(wtDist2[,1],wtDist2[,2])
    # tmpWt[j]<-weighted.mean(dist.tmp,mrcaVec)
  }
  # weighted.distgen2[i]<-weighted.mean(wtDist[,1],wtDist[,2])
  weighted.distgen2[i]<-mean(tmpWt)
  weighted.distgen2.mrcaTRUE[i]<-mean(tmpWt2)
}
sub.modfit[,post]<-weighted.distgen2
true.modfit[,post]<-distgenfromMRCA2
print(post)
}
## Distance#####################
GenSeqMax<-seq(1,50,1)
GenSeqMin<-GenSeqMax-2
GenSeqMin[which(GenSeqMin<0)]<-0
GenSeqMid<-(GenSeqMin+GenSeqMax)/2

distByGenStart2<-distByGenStart1<-distByGenSamp<-distByGen<-rep(NaN,length(GenSeqMax))
for(i in 1:length(GenSeqMin)){
  a<-which(dat.in.all$meanGen>=GenSeqMin[i]&dat.in.all$meanGen<GenSeqMax[i])
  distByGen[i]<-mean(dat.in.all$dist[a])
  
  a<-which(dat.in.all$nogens>=GenSeqMin[i]&dat.in.all$nogens<GenSeqMax[i])
  distByGenStart1[i]<-mean(dat.in.all$distStart1[a])
  
  a<-which(dat.in.all$nogens2>=GenSeqMin[i]&dat.in.all$nogens2<GenSeqMax[i])
  distByGenStart2[i]<-mean(dat.in.all$distStart2[a])
  
  a<-which(dat.in.tmp$meanGen>=GenSeqMin[i]&dat.in.tmp$meanGen<GenSeqMax[i])
  distByGenSamp[i]<-mean(dat.in.tmp$dist[a])
}

data.sim <- cbind(GenSeqMid,distByGen,distByGenSamp,distByGenStart1,distByGenStart2)
data.sim <- data.table(data.sim)
colnames(data.sim) <- c("GenSeqMid","distByGen","distByGenSamp","distByGenStart1","distByGenStart2")

modfit <- cbind(c(1:maxGen),t(apply(sub.modfit, 1, function(x) quantile(x,probs=c(0.025,0.5,0.975),na.rm = TRUE))), 
      t(apply(true.modfit, 1, function(x) quantile(x,probs=c(0.025,0.5,0.975),na.rm = TRUE))))
# modfit <- cbind(c(1:maxGen),distgenfromMRCA2,weighted.distgen,weighted.distgen2,weighted.distgen2.mrcaTRUE)
# colnames(modfit) <- c("gens","distgenfromMRCA2","weighted.distgen","weighted.distgen2","weighted.distgen2.mrcaTRUE")

modfit <- data.table(modfit)
colnames(modfit) <- c("gens","weighted.distgen2.lower","weighted.distgen2","weighted.distgen2.upper",
                      "distgenfromMRCA2.lower","distgenfromMRCA2","distgenfromMRCA2.upper")

simFit <-ggplot()+
  geom_point(data=data.sim,aes(x=GenSeqMid,y=distByGen,color="black"),color="black",alpha=0.7) +
  geom_point(data=data.sim,aes(x=GenSeqMid,y=distByGenSamp,color="red"),color="red",alpha=0.7)+

  geom_line(data=modfit,aes(x=gens,y=distgenfromMRCA2,color="black"),lwd=1.5,color="black",linetype="dashed")+
  geom_ribbon(data=modfit,aes(ymin=weighted.distgen2.lower,ymax=weighted.distgen2.upper,x=gens),alpha=0.2,fill="red")+
  geom_line(data=modfit,aes(x=gens,y=weighted.distgen2,color="red"),lwd=1.5,color="red",linetype="dashed")+
  geom_ribbon(data=modfit,aes(ymin=distgenfromMRCA2.lower,ymax=distgenfromMRCA2.upper,x=gens),alpha=0.4,fill="grey")+
  
  # geom_line(data=modfit,aes(x=gens,y=weighted.distgen2.mrcaTRUE,color="blue"),lwd=1.5,color="blue",linetype="dashed")+
  # geom_line(data=modfit,aes(x=gens,y=weighted.distgen,color="orange"),lwd=1.5,color="orange",linetype="dashed")+
  
  labs(color="Type")+
  theme_bw()+
  scale_x_continuous(limits = c(0, 50), name = "Evolutionary Time (generations)",    # Features of the first axis
                     sec.axis = sec_axis( trans=~./6, name="Evolutionary Time (years)") )+
  theme(axis.text=element_text(size=20), axis.title = element_text(size=20))+
  ylab("Spatial Distance (km)")

simFit
ggsave(simFit,file="./MCMC_model/Simulations/Plots/simFit.pdf")

