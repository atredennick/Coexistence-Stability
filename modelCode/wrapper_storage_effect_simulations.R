##  Wrapper for running simulations of the storage effect
##    model with variable resource.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

# clear the workspace
rm(list=ls())

####
#### Initial conditions and global variables ------------------------
####
maxTime <- 2500 
burn.in <- 1000
DNR <- c(D=c(1,1),N=c(1,1),R=100)
Rmu <- 2                                #mean resource pulse (on log scale)
Rsd <- seq(0,3,by=0.25)                 #std dev of resource pulses (on log scale)
sigE <- c(0, 0.4, 1, 2.5, 5, 7.5, 10)   #environmental cue variability
rho <- seq(-1,1,by=0.25)                #environmental cue correlation between species

resource_sims <- length(Rsd)
env_cue_sims <- length(sigE)
rho_sims <- length(rho)
sims_per_level <- 100