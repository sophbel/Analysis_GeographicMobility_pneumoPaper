likFunc.sim<-function(par){
  
  if(singPar==TRUE){
  extTranMatDat.tmp$pars$homeSus<-par[1]}else{
  ##### Put this back to fit more
  extTranMatDat.tmp$pars$homeSus<-par[1]
  nInfecLoc<-rep(1,9)
  nInfecLoc[1:8]<-exp(par[2:9])
  nInfecLoc<-nInfecLoc/sum(nInfecLoc)}
  # ################################
  
  ### Create mobility matrix for susceptible individuals
  tmp.pHome<-extTranMatDat.tmp$popbyCell
  tmp.pHome<-tmp.pHome/sum(tmp.pHome)
  
  ### Create mobility matrix for susceptible individuals
  tmpbase<-cdr.mat
  tmppar1 <- exp(extTranMatDat.tmp$pars$homeSus)/(1+exp(extTranMatDat.tmp$pars$homeSus))
  tmppar <- min.range+tmppar1*(max.range-min.range)
  tmpdiag<-diag(tmpbase)-tmppar1
  tmpdiag[which(tmpdiag>0.99999)]<-0.99999
  diag(tmpbase)<-0
  tmpbase<-sweep(tmpbase,1,rowSums(tmpbase),"/")
  tmpbase<-sweep(tmpbase,1,(1-tmpdiag)/(1-diag(tmpbase)),"*")
  diag(tmpbase)<-tmpdiag
  tmp.sick<-tmp<-tmpbase
  
 
  move3<-tcrossprod(tmp.sick,tmp.sick)
  move4<-sweep(move3,2,tmp.pHome,"*")
  TranMat.tmp2<-TranMat.tmp<-sweep(move4,1,rowSums(move4),"/")

  # gensA<-gensB<-(1:maxGen)
  maxgen.tmp<-maxGen
  TranMatArray<-array(NA,c(nlocs,nlocs,maxgen.tmp))
  TranMatArray[,,1]<-TranMat.tmp2
  for (j in 2:maxgen.tmp){TranMatArray[,,j]<-eigenMapMatMult(TranMatArray[,,j-1], TranMat.tmp2)
  }
  
  
  mrcaVec<-extTranMatDat.tmp$popbyCell
  # mrcaVec<-tmp.pHome

  llIndPair<-function(ii){
    
    probByGenA<-probByGen.tmp[ii,1:maxGen,1]
    probByGenB<-probByGen.tmp[ii,1:maxGen,2]
    
    minProb=0
    
    gensA<-(1:(maxGen))[which(probByGenA>minProb)]
    gensB<-(1:(maxGen))[which(probByGenB>minProb)]
    

    ndetect<-hist(dat.in2$loc1,plot=FALSE,breaks=seq(0.5,9.5,1))$counts
    probByloc1<-ndetect/nInfecLoc

    
    ndetect<-hist(dat.in2$loc2,plot=FALSE,breaks=seq(0.5,9.5,1))$counts
    probByloc2<-ndetect/nInfecLoc

    
    TranMatArrayA<-TranMatArray[,,gensA]
    if (length(gensA)>1) { TranMatArrayA<-sweep(TranMatArrayA,3,probByGenA[gensA],"*")} else {TranMatArrayA<-TranMatArrayA*probByGenA[gensA]}
    TranMatArrayA<-sweep(TranMatArrayA,2,probByloc1,"*")
    
    TranMatArrayB<-TranMatArray[,,gensB]
    if (length(gensB)>1) { TranMatArrayB<-sweep(TranMatArrayB,3,probByGenB[gensB],"*")} else {TranMatArrayB<-TranMatArrayB*probByGenB[gensB]}
    TranMatArrayB<-sweep(TranMatArrayB,2,probByloc2,"*")
    TranMatArrayB.2<-sweep(TranMatArrayB,1,mrcaVec,"*")
    
    TranMatArrayA2<-sweep(TranMatArrayA,1,mrcaVec,"*")
    TranMatArrayA3<-matrix(TranMatArrayA2,nlocs,nlocs*length(gensA))
    TranMatArrayB2<-matrix(TranMatArrayB,nlocs*length(gensB),nlocs,byrow=T)
    probAllPrs<-eigenMapMatMult(TranMatArrayB2,TranMatArrayA3)
    den<-sum(probAllPrs)
    
    if (length(gensA)>1) {  TranMatArrayA3<-matrix(TranMatArrayA2[,dat.in2$loc1[ii],],nlocs,length(gensA))} else {TranMatArrayA3<-matrix(TranMatArrayA2[,dat.in2$loc1[ii]],nlocs,length(gensA)) }
    if (length(gensB)>1) {  TranMatArrayB2<-matrix(TranMatArrayB[,dat.in2$loc2[ii],],length(gensB),nlocs,byrow=T)} else { TranMatArrayB2<-matrix(TranMatArrayB[,dat.in2$loc2[ii]],length(gensB),nlocs,byrow=T)}
    
    num<-sum(eigenMapMatMult(TranMatArrayB2,TranMatArrayA3))
    
    lik=log(num)-log(den)
    # if(lik<(-6))lik<-(-6)
    
    if(calcAllProbs){
      outlist<-list(all.probs, lik)
      return(outlist)
    }else{
      return(lik)
    }
  }
  
  n.per.core=ceiling(npairs/ncore)
  ntot=n.per.core*ncore
  
  lines.mat<-matrix(1:ntot,nc=ncore,nr=n.per.core)
  
  func<-function(lines){
    if(calcAllProbs==F){
      out<-rep(NaN,length(lines))
      for (iii in 1:length(lines)){out[iii]<-llIndPair(ii=lines[iii])}
      return(out)}else{
        outprobs<-array(NaN,c(length(lines),nlocs,nlocs))
        out2<-rep(NaN,length(lines))
        for (iii in 1:length(lines)){
          a<-llIndPair(lines[iii])
          outprobs[iii,,]<-a[[1]]
          out2[iii]<-a[[2]]
        }
        outlist<-list(outprobs,out2)
        return(outlist)
      }
  }
  
  if(ncore>1){
    aa<-foreach (jj=1:ncore)%dopar%func(lines=lines.mat[,jj])}else{
      aa<-list(func(lines=lines.mat[,1]))
    }

  aa2<-do.call("c",aa)
  LL<-sum(aa2,na.rm=T)
  LL<-LL
  out=sum(LL)
  # attr(out,"IndLik")<-aa2
  print(out)
  
  return(out)
}
