'%notin%'<-Negate('%in%')
#----- Script to extract Landscan population data -----#
library(maptools)
library(sf)
library(rgdal)
library(raster)
library(sp)
library(ggplot2)
library(rgeos)
library(data.table)
library(tidyr)
library(stringr)
################SET UP RR #############
#--- read in shapefile
Sa_shp <- (file = "./data/shapefiles/gadm36_ZAF_3.shp")
shp <- sf::st_read(Sa_shp)


load("./modelinput_data/pop_municipality.2017LS.RData") 
load("./ModelProjections/data/TranMatArray.234x234.RData")
TranMatArray.1<-TranMatArray.234x234
load("./modelinput_data/cdr.mat.town.one.RData")# # [Mobility_ManyMonthsSA.R] ## Probability of movement between each province normalized to carriage rates for each province 
cdr.mat.town<-cdr.mat.town.one
tn <- rownames(cdr.mat.town)
pop2019.town <-pop2019.town[tn]

tn[grep("emalah",tn)]<-c("emalahleni.EC_Eastern Cape","emalahleni.MP_Mpumalanga")
tn[grep("naled",tn)]<-c("naledi.FS_Free State","naledi.NW_Nort West")
tn1<-tn
for (i in 1:length(tn)){
  tn1[i] <- strsplit(tn,"_")[[i]][1]
}
tn1<-str_to_title(tn1)
data.table(tn1)
# moshaweng is botswana 
shp$NAME_3<-as.character(shp$NAME_3)
shp$NAME_3[grep("Emalah",shp$NAME_3)]<-c("Emalahleni.ec","Emalahleni.mp")
shp$NAME_3[grep("Naled",shp$NAME_3)]<-c("Naledi.fs","Naledi.nw")
shp$NAME_3<-as.factor(shp$NAME_3)
shpnames <- shp$NAME_3
tn1[which(is.na(match(tn1,shpnames)))] <- c("City of Cape Town", "City of Johannesburg", "City of Matlosana", "City of Tshwane","Delmas","Dr JS Moroka","eDumbe","Emnambithi/Ladysmith","Greater Marble Hall",
                                            "Kagisano/Molopo","Kai !Garib","//Khara Hais","!Kheis","KwaDukuza","Maluti a Phofung","Mfolozi","Moshaweng","Port St Johns","Pixley Ka Seme","Tlokwe City Council",
                                            "uMhlathuze","uMlalazi","uMngeni","uMshwathi","UMuziwabantu","UPhongolo")
tn1[which(tn1=="Greater Marble Hall")] <- "Ephraim Mogale"
tn1[which(tn1=="Delmas")] <- "Victor Khanye"
tn1[which(tn1=="Moshaweng")] <- "Joe Morolong"
names(pop2019.town) <- tn1
###################population density calculation

matpop<-data.frame("NAME_3"=names(pop2019.town),"population"=pop2019.town)
shp <- merge(shp,matpop,by="NAME_3")

shp$area<-units::set_units(st_area(shp), km^2)
# shp$area<-st_area(shp)/(1000^2)
shp$pop_density<- shp$population/shp$area

shp.trans <- shp %>%
  st_geometry() %>%
  st_transform(crs = '+proj=aeqd ') 
shp.trans_simpl <- st_simplify(shp.trans, preserveTopology = FALSE, dTolerance = 1000)
# centered at input data for low distortion
shp.tmp <- shp
shp.tmp$geometry <- shp.trans_simpl

###############top population sizes
## GT
j_lalo <- data.table(t(data.table(c(-26.2041, 28.0473))))## Joburg
colnames(j_lalo) <- c("latitude","longitude")
j_lalo <- st_as_sf(j_lalo, coords = c("longitude", "latitude"), 
                   crs = 4326, agr = "constant")

## GT ### Not Capital 25.6051° S, 28.3929° E
tshw_lalo <- data.table(t(data.table(c(-25.6051, 28.3929))))## City of Tshwane
colnames(tshw_lalo) <- c("latitude","longitude")
tshw_lalo <- st_as_sf(tshw_lalo, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

## GT ### Not Capital 26.1777° S, 28.3462° E
ekur_lalo <- data.table(t(data.table(c(-26.1777, 28.3462))))## Ekurhuleni
colnames(ekur_lalo) <- c("latitude","longitude")
ekur_lalo <- st_as_sf(ekur_lalo, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

## CT
ct_lalo <- data.table(t(data.table(c(-33.9249, 18.4241)))) ## Cape Town
colnames(ct_lalo) <- c("latitude","longitude")
ct_lalo <- st_as_sf(ct_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
## NC
sp_lalo <- data.table(t(data.table(c(-28.7553, 24.6668)))) ## Sol Plaatjie
colnames(sp_lalo) <- c("latitude","longitude")
sp_lalo <- st_as_sf(sp_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##NW
# 25.8560° S, 25.6403° E
nw_lalo <- data.table(t(data.table(c(-25.8560,25.6403))))
colnames(nw_lalo) <- c("latitude","longitude")
nw_lalo <- st_as_sf(nw_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##MP
# 25.4753° S, 30.9694° E
mp_lalo <- data.table(t(data.table(c(-25.4753,30.9694))))
colnames(mp_lalo) <- c("latitude","longitude")
mp_lalo <- st_as_sf(mp_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##FS
# 29.0852° S, 26.1596° E
fs_lalo <- data.table(t(data.table(c(-29.0852,26.1596))))
colnames(fs_lalo) <- c("latitude","longitude")
fs_lalo <- st_as_sf(fs_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##LP
# 23.8962° S, 29.4486° E
lp_lalo <- data.table(t(data.table(c(-23.8962,29.4486))))
colnames(lp_lalo) <- c("latitude","longitude")
lp_lalo <- st_as_sf(lp_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")
##EC
# 32.9344° S, 27.6435° E
ec_lalo <- data.table(t(data.table(c(-32.9344,27.6435))))
colnames(ec_lalo) <- c("latitude","longitude")
ec_lalo <- st_as_sf(ec_lalo, coords = c("longitude", "latitude"), 
                    crs = 4326, agr = "constant")


##EC ### Port Elizabeth - not the capital Nelson Mandela Bay
# 33.7452° S, 25.5681° E
nmb_lalo <- data.table(t(data.table(c(-33.7452,25.5681))))
colnames(nmb_lalo) <- c("latitude","longitude")
nmb_lalo <- st_as_sf(nmb_lalo, coords = c("longitude", "latitude"), 
                     crs = 4326, agr = "constant")

##KZN ### Ethekwini ### Durban
# 29.8587° S, 31.0218° E
kzn_lalo <- data.table(t(data.table(c(-29.8587,31.0218))))
colnames(kzn_lalo) <- c("latitude","longitude")
kzn_lalo <- st_as_sf(kzn_lalo, coords = c("longitude", "latitude"), 
                     crs = 4326, agr = "constant")



shapes=c("Capitals"=8,"Population >1million"=18)
###############################################################################################
##################### Risk ratio of being at each municipality at 1 year#########
   
  nboot=600000
   ngens=10
   noPop=TRUE
   df.whereat<-matrix(nrow=1,ncol=nboot)
   for( boot in 2:nboot) {
     
     if(noPop==TRUE){specificStart <- sample(234,1)}else
     {specificStart <- sample(234,1,prob=c(pop2019.town)) }
     # specificStart <- 91
     
     df.whereat[1,1] <- wherenext<- specificStart
     # for (gen in 2:ngens){
       df.whereat[1,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,10])
       
     # }
     print(boot)
   }
   
   vecloc <- vector(mode="numeric",length=234)
   df.whereat.tab<-table(df.whereat)
   for (i in 1:234) { 
      vecloc[i] <- (df.whereat.tab[i]/nboot)/mean(df.whereat.tab/nboot)
     # vecloc[i] <- (table(df.whereat)[i]/nboot)/mean(table(df.whereat)/nboot)
     print(i)
   }
   
   vecloc <- data.table(vecloc)
   vecloc$NAME_3 <- tn1
   colnames(vecloc) <- c("Risk_1year","NAME_3")
   vecloc.pop <- cbind(vecloc,pop2019.town)
   
  
   
   shp_simple <- merge(shp.tmp,vecloc,by="NAME_3")
   
save(j_lalo,file="./ModelProjections/data/j_lalo.RData")
save(tshw_lalo,file="./ModelProjections/data/tshw_lalo.RData")
save(ekur_lalo,file="./ModelProjections/data/ekur_lalo.RData")
save(ct_lalo,file="./ModelProjections/data/ct_lalo.RData")
save(kzn_lalo,file="./ModelProjections/data/kzn_lalo.RData")
save(nmb_lalo,file="./ModelProjections/data/nmb_lalo.RData")

if(noPop==TRUE){shp_simple_noPop<-shp_simple
save(shp_simple_noPop,file="./ModelProjections/data/shp_simple_noPop.RData")}else{save(shp_simple,file="./ModelProjections/data/shp_simple.RData")
}

   ##########plot normal
risk1Year <- ggplot(data=shp_simple)+
     geom_sf(data= shp_simple,aes(fill=log(Risk_1year)), lwd = 0) +
     theme(legend.position = "none")+
  
  ### Population >1 million
  geom_sf(data=j_lalo,size=6,shape=18,color="black")+
  geom_sf(data=tshw_lalo,size=6,shape=18,color="black")+
  geom_sf(data=ekur_lalo,size=6,shape=18,color="black")+
  geom_sf(data=ct_lalo,size=6,shape=18,color="black")+
  geom_sf(data=kzn_lalo,size=6,shape=18,color="black")+
  geom_sf(data=nmb_lalo,size=6,shape=18,color="black")+
  # scale_shape_manual(values = shapes,breaks=c("Capitals","Population >1million"),limits=c("Capitals","Population >1million"))+
  theme_classic() +
  scale_fill_distiller( palette ="RdBu", direction = -1,breaks=c(-5,0,2.99),
                        labels=c(0.01,1,20),limits=c(min( log(shp_simple$Risk_1year)),log(max(shp_simple$Risk_1year))) )+
     theme(axis.text=element_blank(),
           plot.subtitle = element_text(color = "blue"),
           plot.caption = element_text(color = "Gray60"),
           legend.text=element_text(size=18),
           # legend.position = c(.9,.15),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +
  legend.position = c(.12,.86),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +

     guides(fill = guide_colorbar(title = "Relative Risk",title.position = "top",
                                  title.theme = element_text(size = 18,
                                                             colour = "gray70",angle = 0)))

risk1Year

ggsave("./ModelProjections/plots/RR1Year.map.Pop.pdf",width = 10,height = 10)
upoverall<-vecloc.pop[which(vecloc.pop$Risk_1year>1)]
dim(upoverall)[1]
tmp<-mat.numIncRisk[1:18,]
quantile(upoverall$Risk_1year,probs=c(0.025,0.5,0.975))
mean(upoverall$Risk_1year)
 ##################################################################################



 
#######################################################
################RURAL VS URBAN###########
############################################

####sample population density
densities.tmp<-as.numeric(shp.tmp$pop_density)
names(densities.tmp)<-shp.tmp$NAME_3
densities<-densities.tmp[names(pop2019.town)]
save(densities,file="./ModelProjections/data/densities.RData")

# min.dens<-c(0,50,500)
# max.dens<-c(50,500,5000)
min.dens<-c(0,500)
max.dens<-c(50,5000)
riskpop.den<-list()
risk.list<-list()
testIncRisk=100
mat.numIncRisk<-matrix(ncol=2,nrow=testIncRisk)
for (g in 1:testIncRisk){
for(d in 1:length(min.dens)){
nboot=1000000
ngens=10
df.whereat<-matrix(nrow=1,ncol=nboot)
samp.tmp<-which(densities>=min.dens[d] & densities<max.dens[d])
samp<-sample(samp.tmp,5)
for( boot in 2:nboot) {
   specificStart <- sample(samp,1)
   # specificStart <- 75
   
   df.whereat[1,1] <- wherenext<- specificStart
  # for (gen in 2:ngens){
    # df.whereat[1,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,12])
    df.whereat[1,boot] <- wherenext <- sample(234,1,prob=TranMatArray.1[wherenext,,10])
    
    
  # }
  print(boot)
}

vecloc <- vector(mode="numeric",length=234)
tab.df.whereat<-table(df.whereat)
for (i in 1:234) { 
   vecloc[i] <- (tab.df.whereat[i]/nboot)/mean(tab.df.whereat/nboot) ### risk of being in X municipality compared to anywhere else on average
   # vecloc[i] <- (tab.df.whereat[i]/nboot)/mean(tab.df.whereat[19]/nboot) ## risk of being in X municipality compared to type you started in
   
  print(i)
}

vecloc <- data.table(vecloc)
vecloc$NAME_3 <- tn1
colnames(vecloc) <- c("Risk_1year","NAME_3")
vecloc.pop <- cbind(vecloc,pop2019.town)
risk.list[[d]]<-vecloc.pop
shp_simple <- merge(shp.tmp,vecloc,by="NAME_3")

####plot rural vs. urban
# risk1Year <- ggplot(data=shp_simple)+
#    geom_sf(data= shp_simple,aes(fill=log(Risk_1year)), lwd = 0) +
#    theme(legend.position = "none")+
#    
#    ### Population >1 million
#    geom_sf(data=j_lalo,size=6,shape=18,color="black")+
#    geom_sf(data=tshw_lalo,size=6,shape=18,color="black")+
#    geom_sf(data=ekur_lalo,size=6,shape=18,color="black")+
#    geom_sf(data=ct_lalo,size=6,shape=18,color="black")+
#    geom_sf(data=kzn_lalo,size=6,shape=18,color="black")+
#    geom_sf(data=nmb_lalo,size=6,shape=18,color="black")+
#    # scale_shape_manual(values = shapes,breaks=c("Capitals","Population >1million"),limits=c("Capitals","Population >1million"))+
#    theme_light() +
#    scale_fill_distiller( palette ="RdBu", direction = -1,breaks=c(-5,0,2.99),
#                          labels=c(0.01,1,20),limits=c(min( log(shp_simple$Risk_1year)),log(max(shp_simple$Risk_1year))) )+
#    theme(axis.text=element_blank(),
#          plot.subtitle = element_text(color = "blue"),
#          plot.caption = element_text(color = "Gray60"),
#          legend.text=element_text(size=12),
#          # legend.position = c(.9,.15),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +
#          legend.position = c(.12,.86),legend.background = element_rect(color="darkgrey", size=0.5, linetype="solid"))  +
#    
#    guides(fill = guide_colorbar(title = "Relative Risk",title.position = "top",
#                                 title.theme = element_text(size = 15,
#                                                            colour = "gray70",angle = 0)))
# risk1Year
# riskpop.den[[d]]<-risk1Year
print(d)
risk.rural<-risk.list[[1]]
uprural<-risk.rural[which(risk.rural$Risk_1year>1)]
risk.urban<-risk.list[[2]]
upurban<-risk.urban[which(risk.urban$Risk_1year>1)]
}

mat.numIncRisk[g,1]<-dim(uprural)[1]
mat.numIncRisk[g,2]<-dim(upurban)[1]
}
# library(patchwork)
# riskpop.den[[1]]+riskpop.den[[2]]+riskpop.den[[3]]
# riskpop.den[[1]]+riskpop.den[[2]]
tmp<-mat.numIncRisk[1:18,]
quantile(tmp[,1],probs=c(0.025,0.5,0.975))
mean(tmp[,1])
quantile(tmp[,2],probs=c(0.025,0.5,0.975))
mean(tmp[,2])
uprural[which(uprural$NAME_3%in%upurban$NAME_3)]

ggsave("./ModelProjections/plots/ruralurban.risk1year.pdf",width=25,height=7)



###general populations
ggplot(data=shp.tmp)+
   geom_sf(data= shp.tmp,aes(fill=as.numeric(pop_density)), lwd = 0)+
   scale_fill_viridis_c(option = "plasma", trans = "sqrt",
                        breaks=c(50,500,2000),limits=c(min(shp.tmp$pop_density),max(shp.tmp$pop_density))) +
   theme_light() +
   theme(axis.text=element_blank(),
         plot.subtitle = element_text(color = "blue"),
         plot.caption = element_text(color = "Gray60"),
         legend.text=element_text(size=20))  +
   labs(fill="Population Density\n(person/km^2)")


# load("/Users/sb62/Documents/Migration/SA_Migration_110422/ModelProjections/RR1Year_popDensities.RData")
