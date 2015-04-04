##  Semi-discrete relative nonlinearity model of two species 
##    coexisting on one essential resource                

##  Authors:  Andrew Tredennick, Peter Adler, and Fred Adler
##  Email:    atredenn@gmail.com
##  Date:     3.19.2015
##  Update:   4.2.2015 -- Implements looping over odSolve for seasons.
##                     -- Makes resource uptake a Hill function, and a
##                        function that can be called within the model.
##                     -- Includes option for a temporally autocorrelated
##                        resource pulse function.
##                     -- Adds 'Rstart' variable to saved output to track
##                        initial condition of resource for each season.

##  MODEL DESCRIPTION
# Species "split" themselves between a dormant, low-mortaility stage (D) and 
# a higher mortality, high growth stage (N).
# The single resource is R.
# There are two sources of variability: an environmental cue that drives the
# storage effect, and resource variability.

# clear the workspace
rm(list=ls())

####
#### Initial conditions, global variables, and parameters ------------------------
####
temporal_autocorrelation <- F       # turn temporal autocorrelation on(T)/off(F)
seasons <- 5                     # number of seasons to simulate
seasons_to_exclude <- 1           # initial seasons to exclude from plots
days_to_track <- 5                  # number of days to recover from odSolve
DNR <- c(D=c(1,1),N=c(1,1),R=10)    # initial conditions
Rmu <- 2                            # mean resource pulse (on log scale)
Rsd <- 0.5                            # std dev of resource pulses (on log scale)
sigE <- 2                           # environmental cue variability
rho <- 0                            # environmental cue correlation between species
parms <- list(
  r = c(10,10),                     # max growth rate for each species
  alpha = c(10,0.5),                 # rate parameter for Hill function 
  beta = c(50,500),                  # shape parameter for Hill function
  mN = c(0.5,0.5),                  # live biomass loss (mortality) rates 
  mD = c(0.001, 0.001)              # dormant biomass loss (mortality) rates
)

####
#### Load relevant libraries ----------------------------------------
####
library(deSolve); library(mvtnorm)


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
days <- c(1:days_to_track)

# Get "germination" fractions for each year
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}

# Loop over seasons
nms <- names(DNR)
gVec <- getG(sigE = sigE, rho = rho, nTime = seasons)
gVec[] <- 0.5

if(temporal_autocorrelation==FALSE){
  Rvector <- rlnorm(seasons,Rmu,Rsd) 
}

if(temporal_autocorrelation==TRUE){
  if(Rsd != 0){
    Rvector <- filter(rlnorm(seasons,Rmu,Rsd), filter=rep(1,3), circular=TRUE)
  }
  if(Rsd == 0){
    Rvector <- rlnorm(seasons,Rmu,Rsd) 
  } 
}

save_seasons <- data.frame(time=NA,D1=NA,D2=NA,N1=NA,N2=NA,R=NA,Rstart=NA,season=NA)

for(season_now in 1:seasons){
  Rstart <- as.list(DNR)$R
  output <- as.data.frame(ode(y = DNR, times = days, 
                              func = updateDNR, parms = parms))
  DNR <- as.numeric(output[nrow(output),nms])
  names(DNR) <- nms
  R_update <- Rvector[season_now]
  DNR <- update_DNR(season_now, DNR, gVec[season_now,])
  names(DNR) <- nms
  new_season <- as.data.frame(output)
  new_season$season <- season_now
  new_season$Rstart <- Rstart
  save_seasons <- rbind(save_seasons, new_season)
}


####
####  Make some plots ----------------------------------
####
#take out first couple seasons
save_seasons <- subset(save_seasons, season>seasons_to_exclude)
matplot(seq_along(along.with = save_seasons$time), 
        save_seasons[,c("N1","N2")], type="l")

#plot just seasons
end_of_season <- which(save_seasons$time == max(days))
end_seasons <- save_seasons[end_of_season,]
matplot(end_seasons$season, end_seasons[,c("N1","N2")], type="l")
plot(N1 ~ Rstart, end_seasons)
plot(N2 ~ Rstart, end_seasons)
plot(N1 ~ N2, end_seasons)

#Plot the resource uptake function
# R <- seq(0,100,1)
# out_r <- matrix(ncol=2, nrow=length(R))
# alpha <- c(10,10)
# beta <- c(50,50)
# for(i in 1:nrow(out_r)){
#   out_r[i,1] <- uptake_R(5, R[i], alpha[1], beta[1])
#   out_r[i,2] <- uptake_R(5, R[i], alpha[2], beta[2])
# }
# matplot(R, out_r, type="l")

#Plot histograms of resource supply rates at different sigRs
# sigR <- seq(0,1,0.25)
# par(mfrow=c(2,3))
# for(i in 1:length(sigR)){
#   hist(rlnorm(100,2,sigR[i]))
# }
# plot(filter(rlnorm(100,2,1), filter=rep(1,3), circular=TRUE))

R <- seq(0,100,1)
out_r <- matrix(ncol=2, nrow=length(R))
alpha <- c(10,0.5)
beta <- c(50,500)
for(i in 1:length(out_r)){
  out_r[i,1] <- uptake_R(10, R[i], alpha[1], beta[1])
  out_r[i,2] <- uptake_R(10, R[i], alpha[2], beta[2])
}
matplot(R, out_r, type="l")

