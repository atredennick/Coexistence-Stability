##  Semi-discrete storage effect model of two species 
##  coexisting on one essential resource.                

##  Authors:  Andrew Tredennick, Peter Adler, and Fred Adler
##  Email:    atredenn@gmail.com
##  Date:     4.24.2015 -- Copied over from relative_nonlinearity.R and 
##                         variable germination included for storage effect.
##            5.21.2016 -- Updated model to have most dynamics occur at the
##                         end of season, discrete time intervals


############################  MODEL DESCRIPTION  ##############################
# Species "split" themselves between a dormant, low-mortaility stage (D) and  #
# a higher mortality, high growth stage (N).                                  #
# The single resource is R.                                                   #
# There are two sources of variability: an environmental cue that drives the  #
# storage effect, and resource variability.                                   #
###############################################################################

# Clear the Workspace
rm(list=ls())

ptm <- proc.time() # start time


####
#### Initial Conditions, Global Variables, and Parameters ----------------------
####
seasons <- 1000                  # number of seasons to simulate
seasons_to_exclude <- 200        # initial seasons to exclude from plots
days_to_track <- 20              # number of days to simulate in odSolve
DNR <- c(D=c(1,1),N=c(1,1),R=10) # initial conditions
Rmu <- 3                         # mean resource pulse (on log scale)
Rsd_annual <- 0                  # std dev of resource pulses (on log scale)
sigE <- 2                        # environmental cue variance
rho <- 0                         # environmental cue correlation between species

# Within-season parameters
parms <- list(
  r = c(5,5),                    # max growth rate for each species
  alpha = c(5,5),                # rate parameter for Hill function 
  beta = c(20,20),               # shape parameter for Hill function
  eps = c(0.2,0.19)               # resource-to-biomass efficiency
)

# End-of-season transition parameters
a1 <- 0.50                       # live-to-dormant biomass fraction; spp1
a2 <- 0.50                       # live-to-dormant biomass fraction; spp2
b1 <- 0                          # adult survivorship; spp1 (0 if annual, >0 if perennial)
b2 <- 0                          # adult survivorship; spp2 (0 if annual, >0 if perennial)
m1 <- 0.1                        # dormant mortality; spp1
m2 <- 0.1                        # dormant mortality; spp2
theta1 <- 0                      # resource recycling fraction; spp1
theta2 <- 0                      # resource recycling fraction; spp2
c <- 0                           # resource carry-over fraction

####
#### Load relevant libraries ----------------------------------------
####
library('deSolve')
library('mvtnorm')



####
#### Model functions -------------------------------------------------
####
## Continuous model
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = 0
    dD2dt = 0
    dN1dt = N1*eps[1]*(uptake_R(r[1], R, alpha[1], beta[1]))
    dN2dt = N2*eps[2]*(uptake_R(r[2], R, alpha[2], beta[2]))
    dRdt = -1 * (dN1dt + dN2dt)
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt)) #output
  })
}

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, alpha, beta){
  return((r*R^alpha) / (beta^alpha + R^alpha))
}

## Discrete model
update_DNR <- function(t,DNR,gs,a1,a2,m1,m2,b1,b2,theta1,theta2,c){
  with (as.list(DNR),{
    g1 <- gs[1]
    g2 <- gs[2]
    D1 <- a1*N1 + D1*(1-g1)*(1-m1)
    D2 <- a2*N2 + D2*(1-g2)*(1-m2)
    N1 <- b1*(1-a1)*N1 + g1*(D1+(a1*N1))*(1-m1)
    N2 <- b2*(1-a2)*N2 + g2*(D2+(a2*N2))*(1-m2)
    R <- theta1*(1-a1)*N1 + theta2*(1-a2)*N2 + c*R + Rvector[t]
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
Rvector <- rlnorm(seasons, Rmu, Rsd_annual)
saved_outs <- matrix(ncol=5, nrow=seasons)
for(season_now in 1:seasons){
  output <- ode(y = DNR, times=days,
                func = updateDNR, parms = parms)
  DNR <- output[nrow(output),nms]
  saved_outs[season_now,] <- DNR
  names(DNR) <- nms
  DNR <- update_DNR(season_now, DNR, gVec[season_now,],
                    a1=a1,a2=a2,m1=m1,m2=m2,b1=b1,b2=b2,
                    theta1=theta1,theta2=theta2,c=c)
  names(DNR) <- nms
}

matplot(saved_outs[101:seasons,c(3:4)], type="l")


proc.time() - ptm
