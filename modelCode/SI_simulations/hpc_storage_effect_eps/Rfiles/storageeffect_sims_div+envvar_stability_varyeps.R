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
library(deSolve)
library(Matrix)
library(mvtnorm)

##  Get simulation id from SLURM
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
i <- as.numeric(myargument)
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary -- here, initial conditions to vary richness
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))
eps <- rbind(c(eps=c(0.2,0.2,0.2,0.2)),
             c(eps=c(0.8,0.8,0.8,0.8)))

n_sig_e <- 15 # Number of cue variance levels
sig_e_vec <- pretty(seq(0.1, 2, length.out=n_sig_e), n_sig_e) # Make a pretty vector
prm <- expand.grid(as.matrix(sig_e_vec), 1:2, 1:4)

DNR_repped <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
eps_repped <- matrix(rep(eps,each=20),ncol=ncol(eps))

colnames(prm) <- c("sigE", "epsid", "dnr_id")
colnames(DNR_repped) <- colnames(DNR)
colnames(eps_repped) <- colnames(eps)
prm <- cbind(prm, eps_repped, DNR_repped)
prm_full <- subset(prm, select = -c(epsid,dnr_id))


##  Define constant parameters in list
constant_parameters <- list (
  seasons = 5000,                  # number of seasons to simulate
  days_to_track = 100,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 0,                  # std dev of resource pulses (on log scale)
  # sigE = 4,                        # environmental cue variance
  rho = 0,                        # environmental cue correlation between species
  alpha1 = 0.5,                    # live-to-dormant biomass fraction; spp1
  alpha2 = 0.495,                  # live-to-dormant biomass fraction; spp2
  alpha3 = 0.49,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.485,                  # live-to-dormant biomass fraction; spp4
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  eta3 = 0.1,                      # dormant mortality; spp3
  eta4 = 0.1                       # dormant mortality; spp4
)

# Growth function parameters
grow_parameters <- list (
  r = c(0.2,0.2,0.2,0.2),   # max growth rate for each species
  a = c(2,2,2,2),           # rate parameter for Hill function 
  b = c(2.5,2.5,2.5,2.5)   # shape parameter for Hill function
  # eps = c(0.2,0.2,0.2,0.2)  # resource-to-biomass efficiency
)



# Make on long vector of named parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters))


# Add in variable parameters to form parameter matrix
constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(prm_full), 
                                ncol=length(constant_param_vec), byrow = TRUE)
colnames(constant_param_matrix) <- names(constant_param_vec)
parameter_matrix <- cbind(constant_param_matrix, prm_full)


##  Run a single instance
outs <- do.call(simulate_model, as.list(parameter_matrix[i,]))

saveRDS(outs, paste0("storageeffect_div_envvar_stability_hiloweps_", i, ".RDS"))