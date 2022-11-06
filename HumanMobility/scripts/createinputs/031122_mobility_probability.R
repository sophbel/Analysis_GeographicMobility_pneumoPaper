## Facebook Movement Data converted to mobility between provinces for South Africa and normalized to number of facebook users per province----
'%notin%' <- Negate('%in%')
library(stringi)
library(maptools)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

### to load from here
load("./data/facebook/mvment_SA.provinces.RData")
load("./data/landscan2017/landscan_populations.RData")
mvment_SA.provinces.ZA$n_baseline_mean <- mvment_SA.provinces.ZA$n_baseline/17

# ## create mobility matrix Province----
SA_pair_BL.raw <- xtabs(n_baseline~start_province + end_province, mvment_SA.provinces.ZA)
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)

## Normalized by population
pop.prov <- with(landscan_populations, tapply(Population, list(NAME_1), FUN=sum)) ## Sum so can do it for all provinces
SA_pair_BL.tmp <- t(apply(SA_pair_BL.raw.one, 1, function(x) x*(1/(pop.prov[rownames(SA_pair_BL.raw.one)]  ) ) ) ) 
tmp.pop <- rowSums(SA_pair_BL.tmp)
## Forced to one
SA_pair_BL.normal.one <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )
SA_pair_BL.normal.pop <- apply(SA_pair_BL.tmp, 2, function(x) x/tmp.pop  )
probmob.Prov <- SA_pair_BL.normal.one
probmob.Prov.RAW <- SA_pair_BL.tmp
###
# ## create mobility matrix Municipality----
SA_pair_BL.raw <- xtabs(n_baseline~start_location + end_location, mvment_SA.provinces.ZA)##BASELINE
# SA_pair_BL.raw <- xtabs(n_crisis~start_polygon_name + end_polygon_name, mvment_SA.provinces.ZA)##CRISIS
# ##INDEX IDS TO POLYGON NAMES OR CONCATENATE PROVINCE NAMES ON TO POLYGON NAMES
SA_pair_BL.raw.one <-  SA_pair_BL.raw+1
tmp.one <- rowSums(SA_pair_BL.raw.one)
## Normalized by population
landscan_populations$location <- paste0(landscan_populations$NAME_3,"_",landscan_populations$NAME_1)
pop.munic <- landscan_populations$Population
names(pop.munic) <- landscan_populations$location
# SA_pair_BL.tmp <- t(apply(SA_pair_BL.raw.one, 1, function(x) x*(1/(pop.munic[rownames(SA_pair_BL.raw.one)]  ) ) ) ) 
SA_pair_BL.tmp <- t(apply(SA_pair_BL.raw.one, 1, function(x) x*(1/(pop.munic[rownames(SA_pair_BL.raw.one)]  ) ) ) )

tmp.pop <- rowSums(SA_pair_BL.tmp)
## Forced to one
# SA_pair_BL.normal.pop <- apply(SA_pair_BL.tmp, 2, function(x) x/tmp.pop  )
SA_pair_BL.normal.one <- apply(SA_pair_BL.raw.one, 2, function(x) x/tmp.one  )

probmob.Munic <- SA_pair_BL.normal.one
probmob.Munic.RAW <- SA_pair_BL.raw

# ## ## Recreate this for a heatmap with NAs grey
# library(ComplexHeatmap)
# library(circlize)
# ## Create mob.matrix with NAs
# SA_pair_BL.raw.none <- probmob.Munic.RAW
# SA_pair_BL.raw.none[SA_pair_BL.raw.none==0] <- NA
# 
# 
# pop.munic <- landscan_populations$Population
# names(pop.munic) <- landscan_populations$location
# SA_pair_BL.tmp2 <- t(apply(SA_pair_BL.raw.none, 1, function(x) x*(1/(pop.munic[rownames(SA_pair_BL.raw.none)]  ) ) ) )
# probmob.Munic.HM <- SA_pair_BL.tmp2
# ## Build index to order by province
# tab <- unique(mvment_SA.provinces.ZA[,c("start_location","start_province")])
# ## Add an index to each municipality that assigns which province it is part of so I can reorder the heatmap according to province and population
# tab$index <- NA
# provs <- names(table(mvment_SA.provinces.ZA$start_province))
# for (i in provs) {
#   tab$index[which(tab$start_province==i)] <- which( names(sort(pop.prov)) == i) 
# }
# col_fun = colorRamp2(c(40,5, 1), c("#9FB1BC","#5c7584","#192024"))
# munic_order <- tab$start_location[order(tab$index)]
# prv.labels <- names(sort(pop.prov))
# probmob.Munic.HM <- probmob.Munic.HM[munic_order,munic_order]
## Plot Heatmaps ----
# hm.munic <- Heatmap(probmob.Munic.HM, 
#               col =col_fun,
#               row_order = rownames(probmob.Munic.HM), column_order = rownames(probmob.Munic.HM), 
#               name = 'Mobility between provinces',
#               show_column_names = TRUE,
#               row_names_gp = gpar(fontsize=6, col=c(rep("#6290C3",27),rep("#002A32",20),rep("#C2847A",19),rep("#7A3B69",18),rep("#F5DD90",25),rep("#3777FF",39),
#                                                     rep("#FFA69E",24),rep("#C6D8D3",51),rep("#E09F3E",10))),
#               column_names_gp = gpar(fontsize=6,col= c(rep("#6290C3",27),rep("#002A32",20),rep("#C2847A",19),rep("#7A3B69",18),rep("#F5DD90",25),rep("#3777FF",39),
#                                                        rep("#FFA69E",24),rep("#C6D8D3",51),rep("#E09F3E",10))),
#               row_names_side = "left",
#               column_names_side = "top",
#               row_title="Starting municipality",
#               column_title = "Ending municipality",
#               heatmap_legend_param = list(
#                 title = "Mobility between municipalities", at = c(15, 5, 0), 
#                 labels = c("0", "0.5", "1") ))
# ### Heatmap for between provinces
# col_fun = colorRamp2(c(0,15, 1.902640e+01 ), rev(c("#FDE0DD","#FA9FB5","#7A0177")))
# hm.prov <- Heatmap(abs(log(probmob.Prov)), 
#               col =col_fun,
#               row_order = names(sort(pop.prov)), column_order = names(sort(pop.prov)),
#               name = 'Mobility between provinces')


## Save file ----
cdr.mat.town.one <- probmob.Munic
cdr.mat.one <- probmob.Prov
# save(cdr.mat, file ="cdr.mat.RData")
# save(cdr.mat.town,file ="cdr.mat.town.RData")
save(cdr.mat.one, file ="./modelinput_data/cdr.mat.one.RData")
save(cdr.mat.town.one,file ="./modelinput_data/cdr.mat.town.one.RData")
# mvmentSA_surroundingcountries=mvmentSA
# save(mvmentSA_surroundingcountries, file="mvmentSA_surroundingcountries.RData")
pop2019.town<-pop.munic
save(pop2019.town,file="./modelinput_data/pop_municipality.2017LS.RData")
save(hm.prov,file="./data/facebook/mobility.plotProvince.RData")
save(hm.munic,file="./data/facebook/mobility.plotMunic.RData")


