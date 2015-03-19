## Semi-discrete relative nonlinearity model of two species   ## 
## coexisting on one essential resource                       ##

####
#### 3/19/2014
#### atredenn@gmail.com
####

# species "split" themselves between a dormant, low-mortaility stage (D) and 
#   a higher mortality, high growth stage (N)
# The single resource is R
# There are two sources of variability: an environmental cue that drives the
#   storage effect, and resource variability

# clear the workspace
rm(list=ls())

####
#### Load relevant libraries
####
library(deSolve)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(gridExtra)


####
#### Model function
####
#continuous model
#TODO: Add dormant state decay
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dN1dt = N1*(r[1]*exp(-k1[1]*exp(-k2[1]*R)) - mN[1])
    dN2dt = N2*(r[2]*exp(-k1[2]*exp(-k2[2]*R)) - mN[2])
    dRdt = -1* ((dN1dt + mN[1]*N1) + (dN2dt + mN[2]*N2))
    list(c(dN1dt, dN2dt, dRdt)) #output
  })
}

#TODO: make discrete model for D -> N transitions
#TODO: try to implement event function and event data at the same time
# initial conditions
maxTime <- 200 
burn.in <- maxTime/10
DNR <- c(N=c(1,1),R=100)
Rmu <- 2 # mean resource pulse (on log scale)
Rsd <-2 # st dev of resource pulses (on log scale)
####
#### Simulate model
####
simTime <- seq(1,maxTime,by=1)
parms <- list(
  r = c(5,1),  # max growth rate for genotype A and a
  k1 = c(20,10),  # right offset for growth rates 
  k2 = c(0.08,0.4),  # rates at which max is approached
  mN = c(0.5,0.5)   # loss (mortality) rates 
)

Rvector <- rlnorm(maxTime,Rmu,Rsd)
eventdat <- data.frame(var=rep("R",maxTime),time = 1:maxTime,value=Rvector,method=rep("add",maxTime))
#run the model
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                           parms = parms, events = list(data = eventdat)))
head(output)
par(mfrow=c(1,2))
matplot(output[,1],output[,2:3],xlab="Time",ylab="N",type="l")
plot(output[,1],output[,4],xlab="Time",ylab="Resource",type="l")

