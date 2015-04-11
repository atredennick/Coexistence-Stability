##  Semi-discrete relative nonlinearity model of two species 
##    coexisting on one essential resource                

##  Authors:  Andrew Tredennick, Peter Adler, and Fred Adler
##  Email:    atredenn@gmail.com
##  Date:     3.19.2015
##  Update:   4.2.2015  -- Implements looping over odSolve for seasons.
##                      -- Makes resource uptake a Hill function, and a
##                         function that can be called within the model.
##                      -- Includes option for a temporally autocorrelated
##                         resource pulse function.
##                      -- Adds 'Rstart' variable to saved output to track
##                         initial condition of resource for each season.
##  Update:   4.8.2015  -- Adds allocation parameter (a) that continuously
##                         draws some fraction of live biomass (N) and allocates
##                         it to the dormant seed bank.
##  Update:   4.10.2015 -- Resource pulses happen within the season.
##                      -- Tracks average seasonal abundance as response.

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
seasons <- 50                   # number of seasons to simulate
seasons_to_exclude <- 0           # initial seasons to exclude from plots
days_to_track <- 120               # number of days to simulate in odSolve
DNR <- c(D=c(1,1),N=c(1,1),R=10)    # initial conditions
Rmu <- 2                            # mean resource pulse (on log scale)
Rsd <- 2                            # std dev of resource pulses (on log scale)
sigE <- 2                           # environmental cue variability
rho <- 0                            # environmental cue correlation between species
parms <- list(
  r = c(10,10),                     # max growth rate for each species
  alpha = c(10,0.5),                 # rate parameter for Hill function 
  beta = c(50,2500),                  # shape parameter for Hill function
  mN = c(0.5,0.5),                  # live biomass loss (mortality) rates 
  mD = c(0.001, 0.001),              # dormant biomass loss (mortality) rates
  a = c(0.5,0.5)                    # allocation fraction of live biomass to seed bank
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
    dD1dt = N1*a[1] - (mD[1]*D1)
    dD2dt = N2*a[2] - (mD[2]*D2)
    dN1dt = N1*(uptake_R(r[1], R, alpha[1], beta[1]) - mN[1] - a[1])
    dN2dt = N2*(uptake_R(r[2], R, alpha[2], beta[2]) - mN[2] - a[2])
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
    N1 <- D1*g1
    N2 <- D2*g2
    D1 <- D1*(1-g1)
    D2 <- D2*(1-g2)
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
gVec[] <- 1

if(temporal_autocorrelation==FALSE){
  Rvector <- runif(seasons,0.1,4)
  Rvector <- rep(2, seasons)
}

if(temporal_autocorrelation==TRUE){
  if(Rsd != 0){
    Rvector <- filter(rlnorm(seasons,Rmu,Rsd), filter=rep(1,3), circular=TRUE)
  }
  if(Rsd == 0){
    Rvector <- rlnorm(seasons,Rmu,Rsd) 
  } 
}

pulse_events <- matrix(nrow=seasons, ncol=days_to_track)
for(do_season in 1:seasons){
  pulse_events[do_season,] <- rlnorm(days_to_track, Rvector[do_season], Rsd)
}
save_seasons <- data.frame(time=NA,D1=NA,D2=NA,N1=NA,N2=NA,R=NA,Rstart=NA,season=NA)

for(season_now in 1:seasons){
  eventdat <- data.frame(var = rep("R", days_to_track),
                         time = 1:days_to_track,
                         value = pulse_events[season_now,],
                         method = rep("add", days_to_track))
  Rstart <- as.list(DNR)$R
  output <- as.data.frame(ode(y = DNR, times = days, 
                              func = updateDNR, parms = parms,
                              events = list(data = eventdat)))
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
# save_seasons <- subset(save_seasons, season>seasons_to_exclude)
# matplot(seq_along(along.with = save_seasons$time), 
#         save_seasons[,c("N1","N2")], type="l")

#plot just average abundance by season
mean_of_season <- ddply(save_seasons, .(season), summarise,
                       mean_N1 = mean(N1),
                       mean_N2 = mean(N2),
                       mean_D1 = mean(D1),
                       mean_D2 = mean(D2))
par(mfrow=c(1,1),mar=c(4,4,1,1),tcl=-0.2,mgp=c(2,0.5,0))
matplot(mean_of_season$season, mean_of_season[,c("mean_N1","mean_N2")], 
        type="l", xlab="season", ylab="abundance")
matplot(mean_of_season$season, mean_of_season[,c("mean_D1","mean_D2")],
        type="l", xlab="season", ylab="abundance")
# plot(N1 ~ Rstart, end_seasons)
# plot(N2 ~ Rstart, end_seasons)
# plot(N1 ~ N2, end_seasons)

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
# 
# R <- seq(0,100,1)
# out_r <- matrix(ncol=2, nrow=length(R))
# alpha <- c(10,0.5)
# beta <- c(50,2500)
# for(i in 1:nrow(out_r)){
#   out_r[i,1] <- uptake_R(10, R[i], alpha[1], beta[1])
#   out_r[i,2] <- uptake_R(10, R[i], alpha[2], beta[2])
# }
# matplot(R, out_r, type="l")

