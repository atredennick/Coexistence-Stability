##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "The relationship between species richness and ecosystem variability is    ##
##  shaped by the mechanism of coexistence"                                   ##
##                                                                       2017 ##
##============================================================================##

# Last update: 2-22-2017

####
#### STORAGE EFFECT FACTORIAL SIMULATIONS IN MULTISPECIES MODEL
#### ISSUE THIS COMMAND IN SHELL BEFORE STARTING R: export OMP_NUM_THREADS=1 
####

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_storage_effect_function_8species.R"
source(fxnfile)                  # Load the function for the simulations
library(deSolve)
library(mvtnorm)
library(Matrix)

##  Get simulation id from SLURM
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i <- as.numeric(myargument)
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary -- here, initial conditions to vary richness
n <- 8
l <- rep(list(0:1), n)
DNR <- expand.grid(l)
DNR <- DNR[which(rowSums(DNR)!=0),]
DNR <- cbind(DNR,DNR,rep(20,nrow(DNR)))
colnames(DNR) <- c(paste0("D",1:8), paste0("N",1:8), "R")

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 5000,                  # number of seasons to simulate
  days_to_track = 100,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
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
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters))


# Add in variable parameters to form parameter matrix
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(DNR), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
parameter_matrix <- cbind(constant_param_matrix, DNR)


##  Run all parameter combinations in paralell
# Returns list of simulation time series with dims = c(nrow(prm), seasons, length(DNR))
outs <- do.call(simulate_model, as.list(parameter_matrix[i,]))

saveRDS(outs, paste0("storage_effect_8species_local_",i,".RDS"))

