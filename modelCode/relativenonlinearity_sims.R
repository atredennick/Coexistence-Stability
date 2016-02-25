##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

## MAIN SIMULATION FILE FOR CONSUMER-RESOURCE SIMULATIONS

####
#### RELATIVE NONLINEARITY
####

rm(list=ls())                    # Erase the memory
fxnfile <- "relative_nonlinearity_fxns.R"
source(fxnfile)                  # Load the name function for the simulations
require(parallel)                # Load the parallel package, to run the simulations in parallel.
library(deSolve)
library(mvtnorm)

nbcores <- 4 # Set number of cores to match machine
set.seed(123456789)

nrho <- 21
rholist <- rep(1, nrho)
nsd <- 21
rsdlist <- pretty(seq(0, 1, length.out=nsd), nsd)

relnonlin_varvars <- expand.grid(rholist,rsdlist)
names(relnonlin_varvars) <- c("rho", "Rsd")
prm <- unique(relnonlin_varvars)
relnon_outs <- mclapply(seq_len(nrow(prm)), function(i) do.call(simRelNonlinModel, prm[i,]), mc.cores=nbcores)

equil_runs <- sapply(relnon_outs, "[", 1)
invasion_runs <- sapply(relnon_outs, "[", 2)

saveRDS(equil_runs, "../simulationResults/relative_nonlinearity_sims.RDS")
saveRDS(invasion_runs, "../simulationResults/relative_nonlinearity_invasion_sims.RDS")
