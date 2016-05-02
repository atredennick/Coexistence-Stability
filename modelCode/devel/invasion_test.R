##  Test of invasion runs for model

rm(list=ls())                    # Erase the memory
fxnfile <- "../simulate_model_function.R"
source(fxnfile)                  # Load the function for the simulations


##  Define constant parameters in list
constant_parameters <- list (
  seasons = 500,                   # number of seasons to simulate
  days_to_track = 20,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 0.0,                # std dev of resource pulses (on log scale)
  sigE = 3,                        # environmental cue variance
  rho = -0.5,                      # environmental cue correlation between species
  alpha1 = 0.50,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.49,                   # live-to-dormant biomass fraction; spp2
  beta1 = 0,                       # adult survivorship; spp1 (0 if annual, >0 if perennial)
  beta2 = 0,                       # adult survivorship; spp2 (0 if annual, >0 if perennial)
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  theta1 = 0,                      # resource recycling fraction; spp1
  theta2 = 0,                      # resource recycling fraction; spp2
  nu = 0                           # resource carry-over fraction
)

# Growth function parameters
grow_parameters <- list (
  r = c(5,5),                   # max growth rate for each species
  a = c(5,5),                   # rate parameter for Hill function 
  b = c(20,20),                 # shape parameter for Hill function
  eps = c(0.5,0.5)              # resource-to-biomass efficiency
)

# Initial conditions
DNR <- c(D=c(1,1),N=c(1,1),R=20) # initial conditions

# Make on long vector of named parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters), unlist(DNR))

# Run a single simulation
outs <- do.call(simulate_model, as.list(constant_param_vec))

equilibrium_runs <- outs[[1]]
matplot(equilibrium_runs[,c(3,4)], type="l")


invasion_runs <- outs[[2]]
matplot(invasion_runs[,c(1,2)], type="l")
abline(h =1, col="blue")
growth_rate <- log(invasion_runs[,2]) - log(1)
mean(growth_rate)

