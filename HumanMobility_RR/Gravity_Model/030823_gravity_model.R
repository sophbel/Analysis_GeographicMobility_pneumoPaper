#### Gravity model for mobility 
#GM: aij is the attraction force between locations i and j. It is directly proportional to masses (or populations) mi and mj and inversely proportional to the squared distance separating them
# dij. mi*mj*dij^-2=GM. The adjusting parameters are alpha, beta, and gamma, where alpha is the exponent adjusting the from probability, beta is adjusting the to probability, and gamma adjusts the distances
### read in locations
load("./modelinput_data/pairwise_geodist.town.RData" )
load("./modelinput_data/pop_municipality.2017LS.RData")
# n<-50
# xlocs<-runif(n,0,100)
# ylocs<-runif(n,0,100)
# popsize<-runif(n,5000,50000)
popsize<-pop2019.town
alpha=1.2
beta=1.5
gamma=2
# dists<-as.matrix(dist(cbind(xlocs,ylocs)))
dists<-pairwise_geodist.town
diag(dists)<-1
n<-dim(dists)[1]
probMove<-matrix(NaN,n,n)
for(i in 1:n){
  for(j in 1:n){
    probMove[i,j]<-popsize[i]^alpha*popsize[j]^beta/dists[i,j]^gamma
  }
}
# diag(probMove)<-NA
rowtot<-rowSums(probMove)
probMove_grav_munic <- apply(probMove, 2, function(x) x/rowtot  )

timeWindow<-35
probStay<-1-(diag(probMove_grav_munic))^timeWindow
tmp<-probMove_grav_munic
diag(tmp)<-0
tmp<-sweep(tmp,1,rowSums(tmp),"/")
tmp<-sweep(tmp,1,(1-probStay)/(1-diag(tmp)),"*")
diag(tmp)<-probStay
probMoveWeek<-tmp

