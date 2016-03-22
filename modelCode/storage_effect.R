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
  a = c(5,5),                    # rate parameter for Hill function 
  b = c(20,20),                  # shape parameter for Hill function
  eps = c(0.2,0.2)               # resource-to-biomass efficiency
)

# End-of-season transition parameters
alpha1 <- 0.50                   # live-to-dormant biomass fraction; spp1
alpha2 <- 0.49                   # live-to-dormant biomass fraction; spp2
beta1 <- 0                       # adult survivorship; spp1 (0 if annual, >0 if perennial)
beta2 <- 0                       # adult survivorship; spp2 (0 if annual, >0 if perennial)
eta1 <- 0.1                      # dormant mortality; spp1
eta2 <- 0.1                      # dormant mortality; spp2
theta1 <- 0                      # resource recycling fraction; spp1
theta2 <- 0                      # resource recycling fraction; spp2
nu <- 0                          # resource carry-over fraction

####
#### Load relevant libraries ----------------------------------------
####
library('deSolve')
library('mvtnorm')



####
#### Model functions -------------------------------------------------
####
## Continuous model
updateNR <- function(t, NR, parms){
  with(as.list(c(NR, parms)), {
    dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
    dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
    dRdt = -1 * (dN1dt/eps[1] + dN2dt/eps[2])
    list(c(dN1dt, dN2dt, dRdt)) #output
  })
}

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, a, b){
  return((r*R^a) / (b^a + R^a))
}

## Discrete model
update_DNR <- function(t,DNR,gammas,alpha1,alpha2,eta1,eta2,beta1,beta2,theta1,theta2,nu){
  with (as.list(DNR),{
    g1 <- gammas[1]
    g2 <- gammas[2]
    D1 <- alpha1*N1 + D1*(1-g1)*(1-eta1)
    D2 <- alpha2*N2 + D2*(1-g2)*(1-eta2)
    N1 <- beta1*(1-alpha1)*N1 + g1*(D1+(alpha1*N1))*(1-eta1)
    N2 <- beta2*(1-alpha2)*N2 + g2*(D2+(alpha2*N2))*(1-eta2)
    R <- theta1*(1-alpha1)*N1 + theta2*(1-alpha2)*N2 + nu*R + Rvector[t]
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
nmsDNR <- names(DNR)
NR <- DNR[c("N1", "N2", "R")] 
nmsNR <- names(NR)
gVec <- getG(sigE = sigE, rho = rho, nTime = seasons)
Rvector <- rlnorm(seasons, Rmu, Rsd_annual)
saved_outs <- matrix(ncol=5, nrow=seasons)
for(season_now in 1:seasons){
  output <- ode(y = NR, times=days,
                func = updateNR, parms = parms)
  NR <- output[nrow(output),nmsNR]
  DNR <- c(DNR[c("D1","D2")], NR)
  saved_outs[season_now,] <- DNR
  names(DNR) <- nmsDNR
  DNR <- update_DNR(season_now, DNR, gVec[season_now,],
                    alpha1=alpha1,alpha2=alpha2,
                    eta1=eta1,eta2=eta2,
                    beta1=beta1,beta2=beta2,
                    theta1=theta1,theta2=theta2,nu=nu)
  names(DNR) <- nmsDNR
  NR <- DNR[c("N1", "N2", "R")] 
  names(NR) <- nmsNR
}
matplot(saved_outs[101:seasons,c(3:4)], type="l", 
        xlab="Season", ylab="Abundance",
        main=expression(paste(alpha[1],"= 0.5; ", alpha[2],"=0.49; ", sigma[E],"=0.2")))


proc.time() - ptm
