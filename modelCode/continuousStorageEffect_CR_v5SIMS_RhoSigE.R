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
maxTime <- 2000 #simulation run time
burn.in <- maxTime*0.5
nSims <- 5
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
sigEvec <- c(0.1,0.5,1,3,5) #environmental cue variability
rhoVec <- c(1,0.5,0,-0.5,-1) #environmental cue correlation between species
Rmu <- 2
# RsdVec <- c(0,0.25,0.5,0.75,1)
Rsd <- 0

####
#### Storage matrices for output
####
covOut <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
sdSumOut <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
cvOut <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
cvRes <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
meanOut <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
sdOut <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
sppRatio <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))
synchOut <- matrix(nrow=length(sigEvec), ncol=length(rhoVec))

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
    dRdt = a*(S-R) - ((N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R))))) + (N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R))))))
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
  for(sigRep in 1:length(sigEvec)){
  #function for transition fraction time series
  getG <- function(sigE, rho, nTime){
    varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
    e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
    g <- exp(e) / (1+exp(e))
    return(g)
  }
  gVec <- getG(sigE = sigEvec[sigRep], rho = rhoVec[gRep], nTime = maxTime)
  gVec1 <- gVec[,1]
  gVec2 <- gVec[,2]
  
    Rnow <- rlnorm(maxTime, Rmu, Rsd)
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
    covN <- numeric(nSims)
    sdNsum <- numeric(nSims)
    synch <- numeric(nSims)
    
    for(sim in 1:nSims){
      DNR <- c(D=c(50,50), N=c(50,50),R=2) #initial conditions
      #run the model
      output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                                 parms = parms, events = list(func = gfun, times=simTime))) 
#       cvN[sim] <- sd(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])/mean(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      totN <- mean(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      sdNsum[sim] <- sd(output[burn.in:maxTime,4])+sd(output[burn.in:maxTime,5])
      covN[sim] <- cov(output[burn.in:maxTime,4],output[burn.in:maxTime,5])
      sdN[sim] <- sd(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      cvN[sim] <- totN/sdN[sim]
      cvR[sim] <- sd(output[burn.in:maxTime,6])/mean(output[burn.in:maxTime,6]) 
      avgN[sim] <- mean(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])
      synch[sim] <- var(output[burn.in:maxTime,4]+output[burn.in:maxTime,5])/(sd(output[burn.in:maxTime,4])+sd(output[burn.in:maxTime,5]))^2
      ifelse(is.na(sdN[sim])==TRUE, stop("SD is undefined."), "sd is good")
#       rN[sim] <- mean(output[burn.in:maxTime,5])/mean(output[burn.in:maxTime,4])
      rN[sim] <- cov(output[burn.in:maxTime,4],output[burn.in:maxTime,5])
      print(paste("Simulation",sim,"for sigma", sigRep, "with rho", gRep))
    }#end simulation loop
    
    covOut[sigRep,gRep] <- mean(covN)
    sdSumOut[sigRep,gRep] <- mean(sdNsum)
    cvOut[sigRep,gRep] <- mean(cvN)
    cvRes[sigRep,gRep] <- mean(cvR)
    meanOut[sigRep,gRep] <- mean(avgN)
    sdOut[sigRep,gRep] <- mean(sdN)
    sppRatio[sigRep,gRep] <- mean(rN)
    synchOut[sigRep,gRep] <- mean(synch)
    
  }#end environmental cue var reps
}#end rho looping

####
#### Write results to file
####
outList <- list(covOut, sdSumOut, cvOut, cvRes, meanOut, sdOut, sppRatio)
saveRDS(outList, file = "rhoVarysigEvaryResults.rds")

####
#### Plot results
####
colnames(cvRes) <- rhoVec
cvRD <- melt(as.data.frame(cvRes))
colnames(cvOut) <- rhoVec
cvD <- melt(as.data.frame(cvOut))
cvD$sVar <- rep(sigEvec, times=length(rhoVec))
colnames(cvD)[1:2] <- c("eVar", "CV")
cvD$RCV <- cvRD$value
cvD$cvRatio <- cvD$RCV/cvD$CV
colnames(meanOut) <- rhoVec
colnames(sdOut) <- rhoVec
colnames(sppRatio) <- rhoVec
moutD <- melt(as.data.frame(meanOut))
sdoutD <- melt(as.data.frame(sdOut))
ratioD <- melt(as.data.frame(sppRatio))
synchD <- melt(as.data.frame(synchOut))
cvD$meanN <- moutD$value
cvD$sdN <- sdoutD$value
cvD$sppRatio <- ratioD$value
cvD$synch <- synchD$value


g1 <- ggplot(cvD, aes(x=eVar, y=CV, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab(expression(S[T]))+
  scale_color_continuous(name=expression(sigma[E]))+
  theme_bw()

g2 <- ggplot(cvD, aes(x=eVar, y=meanN, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab("Mean Community Biomass")+
  scale_color_continuous(name=expression(sigma[E]))+
  theme_bw()

g3 <- ggplot(cvD, aes(x=eVar, y=sdN, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab("SD of Community Biomass")+
  scale_color_continuous(name=expression(sigma[E]))+
  theme_bw()

# g4 <- ggplot(cvD, aes(x=eVar, y=sppRatio, group=sVar))+
#   geom_line(aes(color=as.numeric(as.character(sVar))))+
#   geom_point(size=8, color="white")+
#   geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
#   xlab(expression(rho))+
#   ylab("cov(N1,N2)")+
#   scale_color_continuous(name=expression(sigma[E]))+
#   theme_bw()
g4 <- ggplot(cvD, aes(x=eVar, y=synch, group=sVar))+
  geom_line(aes(color=as.numeric(as.character(sVar))))+
  geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.numeric(as.character(sVar))))+
  xlab(expression(rho))+
  ylab(expression(phi[C]))+
  scale_color_continuous(name=expression(sigma[E]))+
  theme_bw()

g <- arrangeGrob(g1,g2,g3,g4)
png(filename = "CV_FourPanel_RhoVarySigEVary_StorageEff.png", width = 8, height = 5, units="in", res=200)
print(g)
dev.off()


ggplot(cvD, aes(x=synch, y=log(CV)))+
#   geom_line(aes(color=as.numeric(as.character(sVar))))+
#   geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.character(sVar)))+
  ylab("TS")+
  xlab(expression(phi[C]))+
#   scale_color_continuous(name=expression(sigma[E]))+
  theme_bw()

ggplot(cvD, aes(x=synch, y=log(CV)))+
  #   geom_line(aes(color=as.numeric(as.character(sVar))))+
  #   geom_point(size=8, color="white")+
  geom_point(size=4,aes(color=as.character(eVar)))+
  ylab("TS")+
  xlab(expression(phi[C]))+
  #   scale_color_continuous(name=expression(sigma[E]))+
  theme_bw()
