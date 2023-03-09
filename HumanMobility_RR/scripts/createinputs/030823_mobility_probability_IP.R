## Facebook Movement Data converted to mobility between provinces for South Africa and normalized to number of facebook users per province----
'%notin%' <- Negate('%in%')
library(stringi)
library(maptools)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

### to load from here
load("./data/facebook/raw_data/mvment_SA.provinces.RData")
mvment_SA.provinces.ZA$n_baseline_mean <- mvment_SA.provinces.ZA$n_baseline/17

# ## create mobility matrix Province----
SA_pair_BL.raw <- xtabs(n_baseline~start_province + end_province, mvment_SA.provinces.ZA)
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)

## Forced to one -- THIS IS THE ORGINAL MOBILITY MATRIX
SA_pair_BL.normal.one <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )

### adjusted for infectious period first way
# probMoveMonth_Prov<-1-exp(-SA_pair_BL.normal.one*35)
# tmp.oneIP <- rowSums(probMoveMonth_Prov)
# SA_pair_BL.normal.IP <- apply(probMoveMonth_Prov, 2, function(x) x/tmp.oneIP  )
# probmob.Prov <- SA_pair_BL.normal.IP

### adjusted for infectious period second way
timeWindow<-35
probStay<-1-(diag(SA_pair_BL.normal.one))^timeWindow
tmp<-SA_pair_BL.normal.one
diag(tmp)<-0
tmp2<-sweep(tmp,1,rowSums(tmp),"/")
tmp3<-sweep(tmp2,1,(1-probStay)/(1-diag(tmp2)),"*")
diag(tmp3)<-probStay
probmob.Prov<-tmp3

###
# ## create mobility matrix Municipality----
SA_pair_BL.raw <- xtabs(n_baseline~start_location + end_location, mvment_SA.provinces.ZA)##BASELINE
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)
## Forced to one - THIS IS THE ORGINAL MOBILITY MATRIX
SA_pair_BL.normal.one.munic <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )

### combine probabilities across mean infectious periods
### adjusted for infectious period first way
# probMoveMonth_Munic<-1-exp(-SA_pair_BL.normal.one.munic*35)
# tmp.oneIP <- rowSums(probMoveMonth_Munic)
# SA_pair_BL.normal.IP.munic <- apply(probMoveMonth_Munic, 2, function(x) x/tmp.oneIP  )
# probmob.Munic <- SA_pair_BL.normal.IP.munic

### adjusted for infectious period second way
timeWindow<-35
probStay<-1-(diag(SA_pair_BL.normal.one.munic))^timeWindow
tmp<-SA_pair_BL.normal.one.munic
diag(tmp)<-0
tmp1<-sweep(tmp,1,rowSums(tmp),"/")
tmp2<-sweep(tmp1,1,(1-probStay)/(1-diag(tmp1)),"*")
diag(tmp2)<-probStay
probmob.Munic<-tmp2




## Save file ----
cdr.mat.munic.IP <- probmob.Munic
cdr.mat.IP <- probmob.Prov


save(cdr.mat.IP, file ="./modelinput_data/cdr.mat.IP.RData")
save(cdr.mat.munic.IP,file ="./modelinput_data/cdr.mat.munic.IP.RData")



