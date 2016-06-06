##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

# Last update: 6-06-2016

####
#### STORAGE EFFECT FACTORIAL SIMULATIONS IN MULTISPECIES MODEL
#### ISSUE THIE COMMAND IN SHELL BEFORE STARTING R: export OMP_NUM_THREADS=1 
####

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_model_function_4species.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
n_sig_e <- 11 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
n_rho <- 11 # Number of seasonal standard deviation levels
rho_vec <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vector

##  Create matrix with all possible combinations of varying parameters
varvars <- expand.grid(sig_e_vec, rho_vec )
names(varvars) <- c("sigE", "rho")
prm <- unique(varvars)

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 1000,                  # number of seasons to simulate
  days_to_track = 20,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 0.0,                # std dev of resource pulses (on log scale)
  # sigE = 0,                        # environmental cue variance
  # rho = 1,                         # environmental cue correlation between species
  alpha1 = 0.50,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.49,                   # live-to-dormant biomass fraction; spp2
  alpha3 = 0.48,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.47,                   # live-to-dormant biomass fraction; spp4
  beta1 = 0,                       # adult survivorship; spp1 (0 if annual, >0 if perennial)
  beta2 = 0,                       # adult survivorship; spp2 (0 if annual, >0 if perennial)
  beta3 = 0,                       # adult survivorship; spp3 (0 if annual, >0 if perennial)
  beta4 = 0,                       # adult survivorship; spp4 (0 if annual, >0 if perennial)
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  eta3 = 0.1,                      # dormant mortality; spp3
  eta4 = 0.1,                      # dormant mortality; spp4
  theta1 = 0,                      # resource recycling fraction; spp1
  theta2 = 0,                      # resource recycling fraction; spp2
  theta3 = 0,                      # resource recycling fraction; spp3
  theta4 = 0,                      # resource recycling fraction; spp4
  nu = 0                           # resource carry-over fraction
)

# Growth function parameters
grow_parameters <- list (
  r = c(1,1,1,1),           # max growth rate for each species
  a = c(2,2,2,2),           # rate parameter for Hill function 
  b = c(2.5,2.5,2.5,2.5),   # shape parameter for Hill function
  eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
)

# Initial conditions
DNR <- c(D=c(1,1,1,1),N=c(1,1,1,1),R=20) # initial conditions

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
                   do.call(simulate_model, parameter_matrix[i,])
                 }, 
                 mc.cores=nbcores) # end apply function

saveRDS(outs, "../simulationResults/storage_effect_4species.RDS")

