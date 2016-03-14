##  Semi-discrete relative nonlinearity model of two species 
##    coexisting on one essential resource                

##  Authors:  Andrew Tredennick, Peter Adler, and Fred Adler
##  Email:    atredenn@gmail.com
##  Date:     4.24.2015 -- Copied over from relative_nonlinearity.R and 
##                         variable germination included for storage effect.


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
seasons <- 100                   # number of seasons to simulate
seasons_to_exclude <- 10           # initial seasons to exclude from plots
days_to_track <- 60               # number of days to simulate in odSolve
DNR <- c(D=c(1,1),N=c(1,1),R=10)    # initial conditions
Rmu <- 4                           # mean resource pulse (on log scale)
Rsd <- 0                   # std dev of annual mean resource level
Rsd_annual <- Rsd                            # std dev of resource pulses (on log scale)
sigE <- 4                          # environmental cue variability
rho <- 0                           # environmental cue correlation between species
parms <- list(
  r = c(5,4.9),                     # max growth rate for each species
  alpha = c(5,5),                 # rate parameter for Hill function 
  beta = c(20,20),                  # shape parameter for Hill function
  mN = c(0.5,0.5),                  # live biomass loss (mortality) rates 
  mD = c(0.01, 0.01),              # dormant biomass loss (mortality) rates
  a = c(0.5,0.5)                    # allocation fraction of live biomass to seed bank
)

####
#### Load relevant libraries ----------------------------------------
####
library('deSolve')
library('mvtnorm')
library('plyr')
library('ggplot2')
library('gridExtra')
library('synchrony')


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
    R <- Rvector[t]
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

if(temporal_autocorrelation==FALSE){
  Rvector <- rlnorm(seasons, Rmu, Rsd_annual)
}

# if(temporal_autocorrelation==TRUE){
#   if(Rsd != 0){
#     Rvector <- filter(rlnorm(seasons,Rmu,Rsd), filter=rep(1,3), circular=TRUE)
#   }
#   if(Rsd == 0){
#     Rvector <- rlnorm(seasons,Rmu,Rsd_annual) 
#   } 
# }

pulse_events <- matrix(nrow=seasons, ncol=days_to_track)
for(do_season in 1:seasons){
  pulse_events[do_season,] <- rlnorm(days_to_track, log(Rvector[do_season]), Rsd)
}
pulse_events[,] <- 0
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
save_seasons <- subset(save_seasons, season>seasons_to_exclude)
# matplot(seq_along(along.with = save_seasons$time), 
#         save_seasons[,c("N1","N2")], type="l")

#plot just average abundance by season
mean_of_season <- ddply(save_seasons, .(season), summarise,
                        mean_N1 = mean(N1),
                        mean_N2 = mean(N2),
                        mean_D1 = mean(D1),
                        mean_D2 = mean(D2))


seasonal_total <- apply(mean_of_season[,c("mean_N1","mean_N2")], 1, sum)
matplot(mean_of_season$season[20:100], mean_of_season[20:100,c("mean_N1","mean_N2")], 
        type="l", xlab="season", ylab="abundance", ylim=c(0, max(seasonal_total,na.rm=TRUE)),
        col=c("violet", "grey35"), lty=c(2,2), lwd=c(2,2))
lines(mean_of_season$season, seasonal_total, lwd=3, lty=1)

asynch <- 1 - as.numeric(community.sync(mean_of_season[,c("mean_N1","mean_N2")])[1])
asynch
mean(seasonal_total)/sd(seasonal_total)
