##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

# Last update: 4-27-2016

####
#### RELATIVE NONLINEARITY FACTORIAL SIMULATIONS
#### ISSUE THIE COMMAND IN SHELL BEFORE STARTING R: export OMP_NUM_THREADS=1
####

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_model_function_2species.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
n_rsd <- 11 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector

##  Create matrix with all possible combinations of varying parameters
prm <- as.matrix(rsd_vec)
colnames(prm) <- "Rsd_annual"

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 1000,                  # number of seasons to simulate
  days_to_track = 100,             # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  # Rsd_annual = 0.5,               # std dev of resource pulses (on log scale)
  sigE = 0,                        # environmental cue variance
  rho = 1,                         # environmental cue correlation between species
  alpha1 = 0.50,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.50,                   # live-to-dormant biomass fraction; spp2
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1                       # dormant mortality; spp2
)

# Growth function parameters
grow_parameters <- list (
  r = c(0.2,1),                  # max growth rate for each species
  a = c(2,5),                    # rate parameter for Hill function 
  b = c(2.5,20),                 # shape parameter for Hill function
  eps = c(0.5,0.5)               # resource-to-biomass efficiency
)

# Initial conditions
DNR = c(D=c(1,1),N=c(1,1),R=20) # initial conditions

# Make on long vector of named parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters), unlist(DNR))

# Add in variable parameters to form parameter matrix
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(prm), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
parameter_matrix <- cbind(constant_param_matrix, prm)


##  Run all parameter combinations in paralell
# Returns list of simulation time series with dims = c(nrow(prm), seasons, length(DNR))
outs <- mclapply(seq_len(nrow(prm)), 
                 function(i) {
                   do.call(simulate_model, as.list(parameter_matrix[i,]))
                 }, 
                 mc.cores=nbcores) # end apply function

equilibrium_runs <- sapply(outs, "[", 1)
invasion_runs <- sapply(outs, "[", 2)

saveRDS(equilibrium_runs, "../../simulationResults/relative_nonlinearity_equilibrium_runs.RDS")
saveRDS(invasion_runs, "../../simulationResults/relative_nonlinearity_invasion_runs.RDS")

