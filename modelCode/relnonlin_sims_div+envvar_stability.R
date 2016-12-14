###########################################################################
##  storageeffect_sims_div+envar_stability.R: runs simulations of the    ##
##  multispecies model at different levels of resource variance and all  ##
##  combinations of species.                                             ##
###########################################################################

####
#### ISSUE THIE COMMAND IN SHELL BEFORE STARTING R: export OMP_NUM_THREADS=1 
####

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_model_function_4species.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine/simulations
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary -- here, initial conditions to vary richness
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))

## Define vectors of parameters to vary
n_rsd <- 25 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0.1, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
prm <- as.data.frame(rsd_vec)
colnames(prm) <- "Rsd_annual"


##  Define constant parameters in list
constant_parameters <- list (
  seasons = 5000,                  # number of seasons to simulate
  days_to_track = 100,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  # Rsd_annual = 1.2,                  # std dev of resource pulses (on log scale)
  sigE = 0,                        # environmental cue variance
  rho = 1,                        # environmental cue correlation between species
  alpha1 = 0.5,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.5,                   # live-to-dormant biomass fraction; spp2
  alpha3 = 0.5,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.5,                   # live-to-dormant biomass fraction; spp4
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
  r = c(1,5,10,25)/5,           # max growth rate for each species
  a = c(2,5,10,25),           # rate parameter for Hill function 
  b = c(2.5,20,30,45),   # shape parameter for Hill function
  eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
)




# Make on long vector of named parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters))


# Add in variable parameters to form parameter matrix
DNR_repped <- do.call("rbind", replicate(nrow(prm), DNR,  simplify = FALSE))
prm_repped <- do.call("rbind", replicate(nrow(DNR), prm,  simplify = FALSE))
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(DNR_repped), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
parameter_matrix <- cbind(constant_param_matrix, DNR_repped, prm_repped)


##  Run all parameter combinations in paralell
# Returns list of simulation time series with dims = c(nrow(prm), seasons, length(DNR))
outs <- mclapply(seq_len(nrow(parameter_matrix)), 
                 function(i) {
                   do.call(simulate_model, as.list(parameter_matrix[i,]))
                 }, 
                 mc.cores=nbcores) # end apply function

saveRDS(outs, "../simulationResults/relnonlin_div+envvar_cv_stable2unstable.RDS")

