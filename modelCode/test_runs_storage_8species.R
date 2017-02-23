#########################################################################
##  test_runs.R: does quick simulations of the four species model to   ##
##  test parameter combinations                                        ##
#########################################################################

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_sorage_effect_function.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine/simulations
# set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
# Initial conditions
DNR <- rbind(c(D=rep(1,8),N=rep(1,8),R=20))

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 1000,                  # number of seasons to simulate
  days_to_track = 100,              # number of days to simulate in odSolve
  Rmu = 4,                         # mean resource pulse (on log scale)
  Rsd_annual = 0,                # std dev of resource pulses (on log scale)
  sigE = 10,                        # environmental cue variance
  rho = (-1/8),                         # environmental cue correlation between species
  alpha1 = 0.5,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.495,                   # live-to-dormant biomass fraction; spp2
  alpha3 = 0.49,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.485,                   # live-to-dormant biomass fraction; spp4
  alpha5 = 0.48,                   # live-to-dormant biomass fraction; spp5
  alpha6 = 0.475,                   # live-to-dormant biomass fraction; spp6
  alpha7 = 0.47,                   # live-to-dormant biomass fraction; spp7
  alpha8 = 0.465,                   # live-to-dormant biomass fraction; spp8
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  eta3 = 0.1,                      # dormant mortality; spp3
  eta4 = 0.1,                      # dormant mortality; spp4
  eta5 = 0.1,                      # dormant mortality; spp5
  eta6 = 0.1,                      # dormant mortality; spp6
  eta7 = 0.1,                      # dormant mortality; spp7
  eta8 = 0.1                      # dormant mortality; spp8
)

# Growth function parameters
grow_parameters <- list (
  r = rep(0.2,8),           # max growth rate for each species
  a = rep(2,8),           # rate parameter for Hill function
  b = rep(2.5,8),   # shape parameter for Hill function
  eps = rep(0.7,8)  # resource-to-biomass efficiency
)


# Make on long vector of named parameters
# constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters), unlist(DNR))
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters))


# # Add in variable parameters to form parameter matrix
# constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(rsd_vec), 
#                                 ncol=length(constant_param_vec), byrow = TRUE)
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(DNR), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
parameter_matrix <- cbind(constant_param_matrix, DNR)


##  Run all parameter combinations in paralell
# Returns list of simulation time series with dims = c(nrow(prm), seasons, length(DNR))
outs <- mclapply(seq_len(nrow(DNR)), 
                 function(i) {
                   do.call(simulate_model, as.list(parameter_matrix[i,]))
                 }, 
                 mc.cores=nbcores) # end apply function



# constant_env <- outs[[1]]
# variable_env <- outs[[2]]
# 
# # Look at the time-series
par(mfrow=c(1,1), mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1), las=1)
matplot(outs[[1]][,9:16], type="l", xlab="Growing Season", ylab="Biomass", lty=1)
# matplot(variable_env[,1:4], type="l", xlab="Growing Season", ylab="Biomass", main="Fluctuating Environment", lty=1)
# 
# sd(rowSums(outs[[1]][,5:8])) / mean(rowSums(outs[[1]][,5:8]))
# tail(outs[[1]])
