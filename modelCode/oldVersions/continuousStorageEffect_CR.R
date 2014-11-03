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
maxTime <- 500
c <- c(.5,.5)
b <- c(0.2, 0.2)
mD <- c(0.0001, 0.0001)
r <- c(1, .99)
K <- c(.1, .1)
mN <- c(0.1, 0.1)
a <- 0.5
S <- 5
Rmu <- 1 # mean resource pulse (on log scale)
Rsd <-.6 # st dev of resource pulses (on log scale)
sigE <- c(0.4)
rho <- c(0)

####
#### Model function
####

updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    g1 <- getG1(t)
    g2 <- getG2(t)
    dD1dt = c[1]*b[1]*N1 - g1*D1 - mD[1]*D1
    dD2dt = c[2]*b[2]*N2 - g2*D2 - mD[2]*D2
    dN1dt = N1*(r[1]*R/(K[1]+R)) + g1*D1 - c[1]*b[1]*N1 - mN[1]*N1
    dN2dt = N2*(r[2]*R/(K[2]+R)) + g2*D2 - c[2]*b[2]*N2 - mN[2]*N2
    dRdt = a*(S-R) - (N1*(r[1]*R/(K[1]+R)) + N2*(r[2]*R/(K[2]+R)))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt))
  })
}
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
getG1 <- approxfun(x = c(1:maxTime), y = gVec1, method="linear", rule=2)
getG2 <- approxfun(x = c(1:maxTime), y = gVec2, method="linear", rule=2)

####
#### Simulate model
####
Rvector <- rlnorm(maxTime,Rmu,Rsd)   # random pulses
eventdat <- data.frame(var=rep("R",maxTime),time = 1:maxTime,value=Rvector,method=rep("add",maxTime))
# eventdat2 <- data.frame(var=rep("R",maxTime),time = 1:maxTime,value=Rvector,method=rep("add",maxTime))
# eventdat <- rbind(eventdat1, eventdat2)
  
simTime <- c(1:maxTime)
parms <- list(
  c = c,
  b = b,
#   g = rep(gVecTmp, length.out = maxTime),
  mD = mD,
  r = r,
  K = K,
  mN = mN, 
  a = a,
  S = S
)
DNR <- c(D=c(1,1), N=c(1,1),R=1)
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                           parms = parms, events = list(data = eventdat)))


####
#### Plot results
####
par(mfrow=c(1,3))
matplot(simTime, output[,4:5], type="l", main="Live Biomass")
matplot(simTime, output[,2:3], type="l", main="Dormant Biomass")

R <- c(0:120)
getR <- function(r, R, K){
  out1 <- r[1]*R/(K[1]+R)
  out2 <- r[2]*R/(K[2]+R)  
  return(cbind(out1, out2))
}
tmp<-getR(r, R, K)
matplot(R, tmp, type="l")

