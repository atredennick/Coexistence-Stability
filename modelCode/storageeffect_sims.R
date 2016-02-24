##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

## MAIN SIMULATION FILE FOR CONSUMER-RESOURCE SIMULATIONS

####
#### STORAGE EFFECT
####

rm(list=ls())                    # Erase the memory
fxnfile <- "semi_discrete_consumer_resource_fxns.R"
source(fxnfile)                  # Load the name function for the simulations
require(parallel)                # Load the parallel package, to run the simulations in parallel.


nbcores <- 4 # Set number of cores to match machine
set.seed(123456789)

nrho <- 11
rholist <- pretty(seq(-1, 1, length.out=nrho), nrho)
nsd <- 11
rsdlist <- pretty(seq(0, 1, length.out=nsd), nsd)

nsims <- 20
storage_effect_varvars <- expand.grid(rholist,rsdlist)
names(storage_effect_varvars) <- c("rho", "Rsd")
prm <- storage_effect_varvars
strg_outs <- mclapply(seq_len(nrow(prm)), function(i) do.call(simStorageModel, prm[i,]), mc.cores=nbcores)
saveRDS(strg_outs, "../simulationResults/storage_effect_sims.RDS")


####
####  RELATIVE NONLINEARITY
####
