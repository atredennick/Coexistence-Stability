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
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
n_sig_e <- 5 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 1, length.out=n_sig_e), n_sig_e)
n_rsd <- 5 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0, 1, length.out=n_rsd), n_rsd)

##  Create matrix with all possible combinations of varying parameters
varvars <- expand.grid(sig_e_vec, rsd_vec )
names(varvars) <- c("sigE", "Rsd")
prm <- unique(varvars)

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 500,                   # number of seasons to simulate
  days_to_track = 20,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  # Rsd_annual = 0.5,               # std dev of resource pulses (on log scale)
  # sigE = 0,                        # environmental cue variance
  rho = 1,                         # environmental cue correlation between species
  
  # End-of-season transition parameters
  alpha1 = 0.50,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.50,                   # live-to-dormant biomass fraction; spp2
  beta1 = 0,                       # adult survivorship; spp1 (0 if annual, >0 if perennial)
  beta2 = 0,                       # adult survivorship; spp2 (0 if annual, >0 if perennial)
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  theta1 = 0,                      # resource recycling fraction; spp1
  theta2 = 0,                      # resource recycling fraction; spp2
  nu = 0                          # resource carry-over fraction
)

# Growth function parameters
grow_parameters <- list (
  r = c(5,1),                    # max growth rate for each species
  a = c(5,2),                    # rate parameter for Hill function 
  b = c(20,2.5),                 # shape parameter for Hill function
  eps = c(0.2,0.2)               # resource-to-biomass efficiency
)

# Initial conditions
DNR = c(D=c(1,1),N=c(1,1),R=10) # initial conditions

# Make on long vector of parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters), unlist(DNR))

# Add in variable parameters to form matrix
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(prm), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
constant_param_matrix <- cbind(constant_param_matrix, prm)


##  Run all parameter combinations in paralell
# Returns list of simulation time series with dims = c(nrow(prm), seasons, length(DNR))
relnon_outs <- mclapply(seq_len(nrow(prm)), function(i) do.call(sim_relative_nonlinearity, constant_param_matrix[i,]), mc.cores=nbcores)

# equil_runs <- sapply(relnon_outs, "[", 1)

saveRDS(equil_runs, "../simulationResults/relative_nonlinearity_sims.RDS")
# saveRDS(invasion_runs, "../simulationResults/relative_nonlinearity_invasion_sims.RDS")
