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

### adjusted for infectious period
probMoveMonth_Prov<-1-exp(-SA_pair_BL.normal.one*30)
tmp.oneIP <- rowSums(probMoveMonth_Prov)
SA_pair_BL.normal.IP <- apply(probMoveMonth_Prov, 2, function(x) x/tmp.oneIP  )
probmob.Prov <- SA_pair_BL.normal.IP

###
# ## create mobility matrix Municipality----
SA_pair_BL.raw <- xtabs(n_baseline~start_location + end_location, mvment_SA.provinces.ZA)##BASELINE
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)
## Forced to one - THIS IS THE ORGINAL MOBILITY MATRIX
SA_pair_BL.normal.one.munic <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )

### combine probabilities across mean infectious periods
### adjusted for infectious period
probMoveMonth_Munic<-1-exp(-SA_pair_BL.normal.one.munic*30)
tmp.oneIP <- rowSums(probMoveMonth_Munic)
SA_pair_BL.normal.IP.munic <- apply(probMoveMonth_Munic, 2, function(x) x/tmp.oneIP  )
probmob.Munic <- SA_pair_BL.normal.IP.munic

## Save file ----
cdr.mat.munic.IP <- probmob.Munic
cdr.mat.IP <- probmob.Prov
save(cdr.mat.IP, file ="./modelinput_data/cdr.mat.IP.RData")
save(cdr.mat.munic.IP,file ="./modelinput_data/cdr.mat.munic.IP.RData")



