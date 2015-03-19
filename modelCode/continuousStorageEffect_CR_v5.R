## Semi-discrete temporal storage effect model of two species   ## 
## coexisting on one essential resource                         ##

####
#### 11/4/2014
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
#### Parameters
####
maxTime <- 100 #simulation run time
c <- c(1,1) #not used for now, could be a "cost" parameter for biomass storage
b <- c(0.5, 0.5) #also not used, could be assimilation efficiency
mD <- c(0.0001, 0.0001) #dormant state continuous death rate
r <- c(5, 5) #live state intrinsic growth rates
K2 <- c(0.08, 0.08) #rate of approaching max growth rate
K <- c(10,10) #offset for growth rate function
mN <- c(0.5, 0.5) #live state continuous death rate
a <- 0.5 #resource turnover rate
S <- 50 #average resource supply rate
sVar <- 5 #resource supply rate variability
sigE <- 1 #environmental cue variability
rho <- 0 #environmental cue correlation between species
Rmu <- 5 #mean resource abundance (on log scale)
Rsd <- 0 #temporal standard deviation of resource abundance (on log scale)

####
#### Model function
####
#continuous model
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = -(mD[1]*D1)
    dD2dt = -(mD[2]*D2)
    dN1dt = N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R)))) - mN[1]
    dN2dt = N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R)))) - mN[2]
    dRdt = -1 * ((N1*(r[1]*exp(-K[1]*(exp(-K2[1]*R))))) + (N2*(r[2]*exp(-K[2]*(exp(-K2[2]*R))))))
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
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}
gVec <- getG(sigE = sigE, rho = rho, nTime = maxTime)
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
DNR <- c(D=c(1,1), N=c(1,1),R=100) #initial conditions

#run the model
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR,
                           parms = parms, events = list(func = gfun, times=simTime)))


####
#### Plot results
####
par(mfrow=c(1,3),mar=c(4,4,1,1))
R <- seq(0,100,1)
getR <- function(r, R, K){
  out1 <- (r[1]*exp(-K[1]*(exp(-K2[1]*R))))
  out2 <- (r[2]*exp(-K[2]*(exp(-K2[2]*R)))) 
  return(cbind(out1, out2))
}
tmp<-getR(r, R, K)
synch <- var(output[,4]+output[,5])/(sd(output[,4])+sd(output[,5]))^2
synch

growth_rate <- as.data.frame(rbind(cbind(R, tmp[,1]), cbind(R, tmp[,2])))
growth_rate$species <- c(rep("A", length(R)), rep("B", length(R)))

out_mod <- data.frame(time = rep(simTime, 3),
                      species = c(rep("A", maxTime), rep("B", maxTime), rep("Comm.", maxTime)),
                      biomass = c(output[,4], output[,5], output[,4]+output[,5]))

ts_df <- data.frame(species=c("A", "B", "Comm."),
                    stability=ts)

my_colors <- c("red", "blue")
g1 <- ggplot(growth_rate, aes(x=R, y=V2, color=species))+
  geom_line(size=1)+
  geom_hline(aes(yintercept=mN[1]), linetype=2, color=my_colors[1])+
  geom_hline(aes(yintercept=mN[2]), linetype=2, color=my_colors[2])+
  scale_color_manual(values=my_colors)
g2 <- ggplot(subset(out_mod, species!="Comm."), aes(x=time, y=biomass, color=species))+
  geom_line()+
  scale_color_manual(values=my_colors)
# g3 <- ggplot(ts_df, aes(x=species, y=stability, fill=species))+
#   geom_bar(stat = "identity", width=0.75)
g_all <- grid.arrange(g1, g2, ncol=1)



