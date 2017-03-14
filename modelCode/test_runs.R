#########################################################################
##  test_runs.R: does quick simulations of the four species model to   ##
##  test parameter combinations                                        ##
#########################################################################

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_model_function_4species.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine/simulations
# set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
# Initial conditions
# DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
#              c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
#              c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
#              c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
#              c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20))
set.seed(12345678)

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 5000,                  # number of seasons to simulate
  days_to_track = 100,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 1.2,                # std dev of resource pulses (on log scale)
  sigE = 0,                        # environmental cue variance
  rho = 1,                         # environmental cue correlation between species
  alpha1 = 0.5,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.5,                   # live-to-dormant biomass fraction; spp2
  alpha3 = 0.5,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.5,                   # live-to-dormant biomass fraction; spp4
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  eta3 = 0.1,                      # dormant mortality; spp3
  eta4 = 0.1                      # dormant mortality; spp4
)

# Growth function parameters
# grow_parameters <- list (
#   r = c(1,1,1,1),           # max growth rate for each species
#   a = c(2,2,2,2),           # rate parameter for Hill function
#   b = c(2.5,2.5,2.5,2.5),   # shape parameter for Hill function
#   eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
# )
grow_parameters <- list (
  r   = c(0.2,1.0,2.0,5.0), # max growth rate for each species
  a   = c(2,5,10,25),       # rate parameter for Hill function
  b   = c(2.5,20,30,45),    # shape parameter for Hill function
  eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
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
matplot(outs[[1]][,5:8], type="l", xlab="Growing Season", ylab="Biomass", lty=1)
# matplot(variable_env[,1:4], type="l", xlab="Growing Season", ylab="Biomass", main="Fluctuating Environment", lty=1)

sd(rowSums(outs[[1]][,5:8])) / mean(rowSums(outs[[1]][,5:8]))
tail(outs[[1]])
colMeans(outs[[1]][250:nrow(outs[[1]]),])
