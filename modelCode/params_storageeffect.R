##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

## PARAMETERS FILE, STORAGE EFFECT SIMULATIONS

seasons <- 200                   # number of seasons to simulate
seasons_to_exclude <- 50           # initial seasons to exclude from plots
days_to_track <- 60               # number of days to simulate in odSolve
DNR <- c(D=c(1,1),N=c(1,1),R=10)    # initial conditions
Rmu <- 2                            # mean resource pulse (on log scale)
sigE <- 2                           # environmental cue variability
parms <- list(
  r = c(5,5),                     # max growth rate for each species
  alpha = c(5,5),                 # rate parameter for Hill function 
  beta = c(20,20),                  # shape parameter for Hill function
  mN = c(0.5,0.5),                  # live biomass loss (mortality) rates 
  mD = c(0.01, 0.01),              # dormant biomass loss (mortality) rates
  a = c(0.5,0.5)                    # allocation fraction of live biomass to seed bank
)

nrho <- 11
rholist <- pretty(seq(-1, 1, length.out=nrho), nrho)
nsd <- 11
rsdlist <- pretty(seq(0, 1, length.out=nsd), nsd)

nsims <- 20
storage_effect_varvars <- expand.grid(rholist,rsdlist,c(1:nsims))
names(storage_effect_varvars) <- c("rho", "rsd", "sim")

