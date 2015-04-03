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
maxTime <- 500 
burn.in <- maxTime/10
DNR <- c(D=c(1,1),N=c(1,1),R=100)
Rmu <- 2      #mean resource pulse (on log scale)
Rsd <- 0    #std dev of resource pulses (on log scale)
sigE <- 2     #environmental cue variability
rho <- 0     #environmental cue correlation between species

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
    dN1dt = N1*(uptake_R(r[1], R, alpha[1], beta[1]) - mN[1])
    dN2dt = N2*(uptake_R(r[2], R, alpha[2], beta[2]) - mN[2])
    dRdt = -1 * ((dN1dt + mN[1]*N1) + (dN2dt + mN[2]*N2))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt)) #output
  })
}

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, alpha, beta){
  return((r*R^alpha) / (beta^alpha + R^alpha))
}

R <- seq(0,100,1)
out_r <- matrix(ncol=2, nrow=length(R))
alpha <- c(10,10)
beta <- c(50,50)
for(i in 1:nrow(out_r)){
  out_r[i,1] <- uptake_R(5, R[i], alpha[1], beta[1])
  out_r[i,2] <- uptake_R(5, R[i], alpha[2], beta[2])
}
matplot(R, out_r, type="l")

## Discrete model
update_DNR <- function(t,DNR,gs){
  with (as.list(DNR),{
    g1 <- gs[1]
    g2 <- gs[2]
    D1 <- D1 - (N1+D1)*g1 + N1
    D2 <- D2 - (N2+D2)*g2 + N2
    N1 <- 0+(N1+D1)*g1
    N2 <- 0+(N2+D2)*g2
    R <- R + Rvector[t]
    return(c(D1, D2, N1, N2, R))
  })
}


####
#### Simulate model -----------------------------------------------------
####
simTime <- seq(1,maxTime,by=1)
parms <- list(
  r = c(10,10),          #max growth rate for genotype A and a
  alpha = c(10,10),       #right offset for growth rates 
  beta = c(50,50),    #rates at which max is approached
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


# Set random resource fluctuations


# Run the model
# Loop over seasons
seasons <- 500
days <- c(1:5)
nms <- names(DNR)
gVec <- getG(sigE = sigE, rho = rho, nTime = seasons)
Rvector <- rlnorm(seasons,Rmu,Rsd)
save_seasons <- data.frame(time=NA,D1=NA,D2=NA,N1=NA,N2=NA,R=NA,season=NA)
for(season_now in 1:seasons){
  output <- as.data.frame(ode(y = DNR, times = days, 
                              func = updateDNR, parms = parms))
  DNR <- as.numeric(output[nrow(output),nms])
  names(DNR) <- nms
  R_update <- Rvector[season_now]
  DNR <- update_DNR(season_now, DNR, gVec[season_now,])
  names(DNR) <- nms
  new_season <- as.data.frame(output)
  new_season$season <- season_now
  save_seasons <- rbind(save_seasons, new_season)
}

#take out first couple seasons
save_seasons <- subset(save_seasons, season>100)
matplot(seq_along(along.with = save_seasons$time), save_seasons[,c("N1","N2")], type="l")

#plot just seasons
end_of_season <- which(save_seasons$time == max(days))
end_seasons <- save_seasons[end_of_season,]
matplot(end_seasons$season, end_seasons[,c("N1","N2")], type="l")
# plot(N1 ~ R, end_seasons)
# plot(N2 ~ R, end_seasons)
# plot(N1 ~ N2, end_seasons)


