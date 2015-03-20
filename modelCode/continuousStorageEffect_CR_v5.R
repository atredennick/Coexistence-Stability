## Semi-discrete storage effect model of two species   ## 
## coexisting on one essential resource                ##

####
#### 3/19/2014
#### atredenn@gmail.com
####

##  MODEL DESCRIPTION
# species "split" themselves between a dormant, low-mortaility stage (D) and 
#   a higher mortality, high growth stage (N)
# The single resource is R
# There are two sources of variability: an environmental cue that drives the
#   storage effect, and resource variability

# clear the workspace
rm(list=ls())

####
#### Initial conditions and global variables ------------------------
####
maxTime <- 200 
burn.in <- maxTime/10
DNR <- c(D=c(1,1),N=c(1,1),R=100)
Rmu <- 2      #mean resource pulse (on log scale)
Rsd <- 0    #std dev of resource pulses (on log scale)
sigE <- 1     #environmental cue variability
rho <- 0      #environmental cue correlation between species

####
#### Load relevant libraries ----------------------------------------
####
library(deSolve)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(gridExtra)


####
#### Model function -------------------------------------------------
####
## Continuous model
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = -(mD[1]*D1)
    dD2dt = -(mD[2]*D2)
    dN1dt = N1*(r[1]*exp(-k1[1]*exp(-k2[1]*R)) - mN[1])
    dN2dt = N2*(r[2]*exp(-k1[2]*exp(-k2[2]*R)) - mN[2])
    dRdt = -1* ((dN1dt + mN[1]*N1) + (dN2dt + mN[2]*N2))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt)) #output
  })
}

## Discrete model
gfun <- function(t, y, parms){
  with (as.list(y),{
    g1 <- gVec1[t]
    g2 <- gVec2[t]
    D1 <- D1 - (N1+D1)*g1 + N1
    D2 <- D2 - (N2+D2)*g2 + N2
    N1 <- 0+(N1+D1)*g1
    N2 <- 0+(N2+D2)*g2
    R <- R + Rvector[t]
    return(c(D1, D2, N1, N2, R))
  })
}


#TODO: try to implement event function and event data at the same time

####
#### Simulate model -----------------------------------------------------
####
simTime <- seq(1,maxTime,by=1)
parms <- list(
  r = c(5,5),          #max growth rate for genotype A and a
  k1 = c(20,20),       #right offset for growth rates 
  k2 = c(0.08,0.08),    #rates at which max is approached
  mN = c(0.5,0.5),      #live biomass loss (mortality) rates 
  mD = c(0.001, 0.001) #dormant biomass loss (mortality) rates
)

# Get "germination" fractions for each year
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}
gVec <- getG(sigE = sigE, rho = rho, nTime = maxTime)
gVec1 <- gVec[,1]
gVec2 <- gVec[,2]

# Set random resource fluctuations
Rvector <- rlnorm(maxTime,Rmu,Rsd)

# Run the model
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR, parms = parms,
                           events = list(func = gfun, times=simTime)))

####
#### Make some plots ----------------------------------------------
####
par(mfrow=c(1,2))
matplot(output[,1],output[,2:3],xlab="Time",ylab="N",type="l")
plot(output[,1],output[,4],xlab="Time",ylab="Resource",type="l")
