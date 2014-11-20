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
library(reshape2)
library(ggplot2)
library(gridExtra)

####
#### Parameters
####
maxTime <- 20 #simulation run time
burn.in <- maxTime*0.5
nSims <- 1
c <- c(1,1) #not used for now, could be a "cost" parameter for biomass storage
b <- c(0.5, 0.5) #also not used, could be assimilation efficiency
mD <- c(0.0001, 0.0001) #dormant state continuous death rate
r <- c(2, 2) #live state intrinsic growth rates
K2 <- c(100, 100) #rate of approaching max growth rate
K <- c(25,25) #offset for growth rate function
mN <- c(0.1, 0.1) #live state continuous death rate
a <- 0.5 #resource turnover rate
S <- 10 #average resource supply rate
sVar <- 5 #resource supply rate variability
sigE <- 0.1 #environmental cue variability
rhoVec <- c(1,0.5,0,-0.5,-1) #environmental cue correlation between species
Rmu <- 2
RsdVec <- c(0,0.25,0.5,0.75,1)

####
#### Storage matrices for output
####
cvOut <- matrix(nrow=length(RsdVec), ncol=length(rhoVec))
cvRes <- matrix(nrow=length(RsdVec), ncol=length(rhoVec))
meanOut <- matrix(nrow=length(RsdVec), ncol=length(rhoVec))
sdOut <- matrix(nrow=length(RsdVec), ncol=length(rhoVec))
sppRatio <- matrix(nrow=length(RsdVec), ncol=length(rhoVec))

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
for(gRep in 1:length(rhoVec)){
  #function for transition fraction time series
  getG <- function(sigE, rho, nTime){
    rho <- rhoVec[gRep]
    varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
    e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
    g <- exp(e) / (1+exp(e))
    return(g)
  }
  gVec <- getG(sigE = sigE, rho = rho, nTime = maxTime)
  gVec1 <- gVec[,1]
  gVec2 <- gVec[,2]
  
  for(ResRep in 1:length(RsdVec)){
    Rnow <- rlnorm(maxTime, Rmu, RsdVec[ResRep])
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
    
    #individual sim storage
    cvR <- numeric(nSims)
    cvN <- numeric(nSims)
    avgN <- numeric(nSims)
    sdN <- numeric(nSims)
    rN <- numeric(nSims)
    
    for(sim in 1:nSims){
      DNR <- c(D=c(1500,1500), N=c(1500,1500),R=2) #initial conditions
      #run the model
      output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                                 parms = parms, events = list(func = gfun, times=simTime))) 
      cvN[sim] <- sd(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])/mean(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      cvR[sim] <- sd(output[burn.in:maxTime,6])/mean(output[burn.in:maxTime,6]) 
      avgN[sim] <- mean(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      sdN[sim] <- sd(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      rN[sim] <- mean(output[burn.in:maxTime,5])/mean(output[burn.in:maxTime,4])
      print(c(sim, ResRep, gRep))
    }#end simulation loop
    
    cvOut[ResRep,gRep] <- mean(cvN)
    cvRes[ResRep,gRep] <- mean(cvR)
    meanOut[ResRep,gRep] <- mean(avgN)
    sdOut[ResRep,gRep] <- mean(sdN)
    sppRatio[ResRep,gRep] <- mean(rN)
    
  }#end resource var reps
}#end rho looping


####
#### Plot results
####
colnames(cvRes) <- rhoVec
cvRD <- melt(as.data.frame(cvRes))
colnames(cvOut) <- rhoVec
cvD <- melt(as.data.frame(cvOut))
cvD$sVar <- rep(RsdVec, times=length(rhoVec))
colnames(cvD)[1:2] <- c("eVar", "CV")
cvD$RCV <- cvRD$value
cvD$cvRatio <- cvD$RCV/cvD$CV
colnames(meanOut) <- rhoVec
colnames(sdOut) <- rhoVec
colnames(sppRatio) <- rhoVec
moutD <- melt(as.data.frame(meanOut))
sdoutD <- melt(as.data.frame(sdOut))
ratioD <- melt(as.data.frame(sppRatio))
cvD$meanN <- moutD$value
cvD$sdN <- sdoutD$value
cvD$sppRatio <- ratioD$value


g1 <- ggplot(cvD, aes(x=eVar, y=CV, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab(expression(CV[N[1]+N[2]]))+
  scale_color_continuous(name=expression(sigma[R]))+
  theme_bw()

g2 <- ggplot(cvD, aes(x=eVar, y=meanN, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab("Mean Community Biomass")+
  scale_color_continuous(name=expression(sigma[R]))+
  theme_bw()

g3 <- ggplot(cvD, aes(x=eVar, y=sdN, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab("SD of Community Biomass")+
  scale_color_continuous(name=expression(sigma[R]))+
  theme_bw()

g4 <- ggplot(cvD, aes(x=eVar, y=sppRatio, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab("N2/N1")+
  scale_color_continuous(name=expression(sigma[R]))+
  theme_bw()

g <- arrangeGrob(g1,g2,g3,g4)
png(filename = "CV_FourPanel_RhoVary_StorageEff.png", width = 8, height = 5, units="in", res=200)
g
dev.off()