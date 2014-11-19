#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Semi-discrete temporal storage effect model of two species   ## 
## coexisting on one essential resource                         ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

####
#### 11/4/2014
#### atredenn@gmail.com
####

# species "split" themselves between a dormant, low-mortaility stage (D) and a higher mortality, high growth stage (N)
# the single resource is R
# There are two sources of variability: an environmental cue that drives the storage effect, and resource variability

# clear the workspace
rm(list=ls())

####
#### Load relevant libraries
####
library(deSolve)
library(mvtnorm)

####
#### Parameters
####
maxTime <- 1000 #simulation run time
c <- c(1,1) #not used for now, could be a "cost" parameter for biomass storage
b <- c(0.5, 0.5) #also not used, could be assimilation efficiency
mD <- c(0.0001, 0.0001) #dormant state continuous death rate
r <- c(2.1, 0.15) #live state intrinsic growth rates
K2 <- c(100, 200) #rate of approaching max growth rate
K <- c(25,2) #offset for growth rate function
mN <- c(0.1, 0.1) #live state continuous death rate
a <- 0.5 #resource turnover rate
S <- 5 #average resource supply rate
sVar <- 10 #resource supply rate variability
sigE <- 0 #environmental cue variability
rho <- 0 #environmental cue correlation between species
Rmu <- 2
Rsd=1.5
####
#### Model function
####
#continuous model
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = -(mD[1]*D1)
    dD2dt = -(mD[2]*D2)
    dN1dt = N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R)))) - mN[1]*N1
    dN2dt = N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R)))) - mN[2]*N2
    dRdt = a*(S-R) - (N1*(r[1]*R/(K[1]+R)) + (N2*(r[2]*R/(K[2]+R))))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt)) #output
  })
}

#discrete model: pulses of N-D transitions (difference equations)
gfun <- function(t, y, parms){
  with (as.list(y),{
    g1 <- gVec1[t]
    g2 <- gVec2[t]
    D1 <- D1 - (N1+D1)*g1 + N1
    D2 <- D2 - (N2+D2)*g2 + N2
    N1 <- 0+(N1+D1)*g1
    N2 <- 0+(N2+D2)*g2
    R <- Rnow[t]
    return(c(D1, D2, N1, N2, R))
  })
}


####
#### Simulate model
####
#function for transition fraction time series
gVec1 <- rep(0.2,maxTime)
gVec2 <- rep(0.2,maxTime)

Rnow <- rlnorm(maxTime, Rmu, Rsd)
hist(Rnow)

#forcing for resource supply rate variations
# forcedat <- data.frame(Time = c(1:maxTime),
#                        Supply = rnorm(maxTime,0,sVar))
# forcedat[forcedat$Supply<(-S),2] <- -S #don't let supply go negative
# rand <- approxfun(x = forcedat[,1], y = forcedat[,2], method = "linear", rule = 2)

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
DNR <- c(D=c(0,0), N=c(3000,3000),R=1) #initial conditions

#run the model
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                           parms = parms, events = list(func = gfun, times=simTime)))


####
#### Plot results
####
par(mfrow=c(1,3),mar=c(4,4,1,1))
# matplot(simTime, output[,2:3], type="l", main="Dormant Biomass", 
#         col=c("darkorange", "purple"), xlab="Years (T)", ylab="Biomass (D)")
# plot(simTime, output$R, type="l", main=mean(output$R))
cvs <- c((sd(output[,4])/mean(output[,4])), 
         (sd(output[,5])/mean(output[,5])), 
         (sd(output[,4]+output[,5])/mean(output[,4]+output[,5])))

# plot the growth rate functions for each species
# plot the growth rate functions for each species
R <- seq(0,0.1,0.001)
getR <- function(r, R, K){
  out1 <- (r[1]*exp(-K[1]*(exp(-K2[1]*R))))
  out2 <- (r[2]*exp(-K[2]*(exp(-K2[2]*R)))) 
  return(cbind(out1, out2))
}
tmp<-getR(r, R, K)
matplot(R, tmp, type="l", col=c("darkorange", "purple"), lty=c(1,1), 
        ylab="Instantaneous growth rate (r)", xlab="Resource density (R)")
abline(h = mN[1], lty=2)
matplot(simTime, output[,4:5], type="l", main="Live Biomass", lty=c(1,1), 
        col=c("darkorange", "purple"), xlab="Years (T)", ylab="Biomass (N)")
# matplot(simTime[1900:2000], output[1900:2000,4:5], type="l", main="Live Biomass", lty=c(1,1), 
#         col=c("darkorange", "purple"), xlab="Years (T)", ylab="Biomass (N)")

barplot(cvs, names.arg = c("Spp 1", "Spp 2", "Community"), ylab="C.V.", col=c("darkorange", "purple", "grey"))
box()



