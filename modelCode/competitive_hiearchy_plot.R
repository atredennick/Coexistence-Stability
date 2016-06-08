#########################################################################
##  competitive_hierarchy_plots.R: does two quick simulations of the   ##
##  four species storage effect model to show competitive hierarchies  ##
##  with and without environmental variation                           ##
#########################################################################

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_model_function_4species.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 2 # Set number of cores to match machine/simulations
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
sig_e_vec <- as.data.frame(matrix(c(0,2),2,1))
names(sig_e_vec) <- "sigE"

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 1000,                  # number of seasons to simulate
  days_to_track = 20,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 0.0,                # std dev of resource pulses (on log scale)
  # sigE = 0,                        # environmental cue variance
  rho = -1,                        # environmental cue correlation between species
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
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(sig_e_vec), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
parameter_matrix <- cbind(constant_param_matrix, sig_e_vec)


##  Run all parameter combinations in paralell
# Returns list of simulation time series with dims = c(nrow(prm), seasons, length(DNR))
outs <- mclapply(seq_len(nrow(sig_e_vec)), 
                 function(i) {
                   do.call(simulate_model, parameter_matrix[i,])
                 }, 
                 mc.cores=nbcores) # end apply function

constant_env <- outs[[1]]
variable_env <- outs[[2]]

png("../manuscript/components/fourspecies_comp_hierarchy.png", width = 6, height = 5, units = "in", res = 72)
par(mfrow=c(2,1), mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1), las=1)
matplot(constant_env[,1:4], type="l", xlab="Growing Season", ylab="Biomass", main="Constant Environment", lty=1)
matplot(variable_env[,1:4], type="l", xlab="Growing Season", ylab="Biomass", main="Fluctuating Environment", lty=1)
dev.off()
