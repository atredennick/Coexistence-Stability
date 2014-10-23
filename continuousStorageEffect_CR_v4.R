#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Continuous-time temporal storage effect model of two species ## 
## coexisting on one essential resource                         ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# species "split" themselves between a dormant, low-mortaility stage (D) and a higher mortality, high growth stage (N)
# the single resource is R

rm(list=ls())
library(deSolve)
library(mvtnorm)

####
#### Parameters
####
maxTime <- 1000
c <- c(1,1)
b <- c(0.5, 0.5)
mD <- c(0.0001, 0.0001)
r <- c(2.5, 2.1)
K <- c(.1, .1)
mN <- c(0.1, 0.1)
a <- 0.5
S <- 10
Rmu <- 1 # mean resource pulse (on log scale)
Rsd <-.6 # st dev of resource pulses (on log scale)
sigE <- c(2.5)
rho <- c(-1)

####
#### Model function
####

updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = -(mD[1]*D1)
    dD2dt = -(mD[2]*D2)
    dN1dt = N1*(r[1]*R/(K[1]+R)) - mN[1]*N1
    dN2dt = N2*(r[2]*R/(K[2]+R)) - mN[2]*N2
    dRdt = a*(S-R) - (N1*(r[1]*R/(K[1]+R)) + (N2*(r[2]*R/(K[2]+R))))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt))
  })
}

#function for transition time series
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}
gVec <- getG(sigE = sigE, rho = rho, nTime = maxTime)
gVec1 <- gVec[,1]
gVec2 <- gVec[,2]
# source("getSineTransitions.R")
# gVec1 <- rep(gVecTmp[,1], length.out = maxTime)
# gVec2 <- rep(gVecTmp[,2], length.out = maxTime)


####
#### Simulate model
####
Rvector <- rlnorm(maxTime,Rmu,Rsd)   # random pulses
Rvector <- rep(0,maxTime)
eventdat <- data.frame(var=rep("R",maxTime),time = 1:maxTime,value=Rvector,method=rep("add",maxTime))
# eventdat2 <- data.frame(var=rep("D",maxTime),time = 1:maxTime,value=Rvector,method=rep("mul",maxTime))
# eventdat <- rbind(eventdat1, eventdat2)

gfun <- function(t, y, parms){
  with (as.list(y),{
    g1 <- gVec1[t]
    g2 <- gVec2[t]
    D1 <- D1 - D1*g1 + N1
    D2 <- D2 - D2*g2 + N2
    N1 <- 0+D1*g1
    N2 <- 0+D2*g2
    R <- R
    return(c(D1, D2, N1, N2, R))
  })
}
  
simTime <- seq(1,maxTime,by=1)
parms <- list(
  c = c,
  b = b,
  mD = mD,
  r = r,
  K = K,
  mN = mN, 
  a = a,
  S = S
)
DNR <- c(D=c(30,30), N=c(20,20),R=1)
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                           parms = parms, events = list(func = gfun, times=simTime)))


####
#### Plot results
####
par(mfrow=c(1,3),mar=c(4,4,1,1))
# ylim=c(0, max(output[,4]+output[,5]))
matplot(simTime, output[,4:5], type="l", main="Live Biomass", 
        col=c("darkorange", "purple"), xlab="Years (T)", ylab="Biomass (N)")
# lines(simTime, output[,4]+output[,5], lwd=3, col="dodgerblue")
matplot(simTime, output[,2:3], type="l", main="Dormant Biomass", 
        col=c("darkorange", "purple"), xlab="Years (T)", ylab="Biomass (D)")

R <- seq(0,1,0.01)
getR <- function(r, R, K){
  out1 <- r[1]*R/(K[1]+R)
  out2 <- r[2]*R/(K[2]+R)  
  return(cbind(out1, out2))
}
tmp<-getR(r, R, K)
matplot(R, tmp, type="l", col=c("darkorange", "purple"), 
        ylab="Intrinsic growth rate (r)", xlab="Resource density (R)")

