####
#### Adaptation of Gross 2008 model (just consumer-resource)
####

library(deSolve)

####
#### Parameters
####
r=c(1.5,1.2)
K=c(10,1)
m=c(.05,.05)
NR <- c(N=c(1,1),R=50)
c=c(1,1)
maxDays=120


####
#### Model functions
####
updateNR = function(t, NR, parms) { 
  # Continuous two species resource model with an essential resource
  with(as.list(c(NR, parms)), { 
    dN1dt = ((r[1]*R)/(K[1]+R) - m[1])*N1 ##spp 1 growth
    dN2dt = ((r[2]*R)/(K[2]+R) - m[2])*N2 ##spp 2 growth
    dRdt = -1 * (c[1]*(r[1]*R)/(K[1]+R)*N1 + c[2]*(r[2]*R)/(K[2]+R)*N2) ##resource draw down
    list(c(dN1dt, dN2dt, dRdt))
  }) 
}

getGrowthRate = function(r, K, R){
  out1 <- (r[1]*R)/(K[1]+R)
  out2 <- (r[2]*R)/(K[2]+R)
  return(cbind(out1,out2))
}

####
#### Simulation
####

nSave <- matrix(nrow=maxDays, ncol=2)
days <- 1:maxDays
fNow <- getGrowthRate(r,K,R)
parms <- list(
  r = r,
  K = K,
  m = m
)
output = as.data.frame(ode(y = NR, times = days, func = updateNR, parms = parms))


####
#### Figures
####
par(mfrow=c(1,2))
matplot(days, output[,2:3], type="l")
R <- 0:120
tmp<-getGrowthRate(r=c(2.5,2), K=c(10,1), R=R)
matplot(R,tmp,type="l",ylab="Resource uptake")



