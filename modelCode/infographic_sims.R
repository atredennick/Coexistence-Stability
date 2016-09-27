#########################################################################
##  test_runs.R: does quick simulations of the four species model to   ##
##  test parameter combinations                                        ##
#########################################################################

rm(list=ls())                    # Erase the memory
fxnfile <- "simulate_model_function_4species.R"
source(fxnfile)                  # Load the function for the simulations
require(parallel)                # Load the parallel package

nbcores <- 4 # Set number of cores to match machine/simulations
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
# Initial conditions
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20))
#              c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
#              c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
#              c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
#              c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))

sig_e_vec <- c(0.1, 0.2, 0.75, 1.8) # Make a pretty vector
prm <- as.matrix(sig_e_vec)
colnames(prm) <- "sigE"

##  Define constant parameters in list
constant_parameters <- list (
  seasons = 5000,                  # number of seasons to simulate
  days_to_track = 20,              # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 0.0,                # std dev of resource pulses (on log scale)
  # sigE = 0,                        # environmental cue variance
  rho = 0,                         # environmental cue correlation between species
  alpha1 = 0.50,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.495,                   # live-to-dormant biomass fraction; spp2
  alpha3 = 0.49,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.485,                   # live-to-dormant biomass fraction; spp4
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
# grow_parameters <- list (
#   r = c(1,5,10,25),           # max growth rate for each species
#   a = c(2,5,10,25),           # rate parameter for Hill function 
#   b = c(2.5,20,30,45),   # shape parameter for Hill function
#   eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
# )
grow_parameters <- list (
  r = c(1,1,1,1),           # max growth rate for each species
  a = c(2,2,2,2),           # rate parameter for Hill function 
  b = c(2.5,2.5,2.5,2.5),   # shape parameter for Hill function
  eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
)

# Make on long vector of named parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters), unlist(DNR))
# constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters))


# # Add in variable parameters to form parameter matrix
# constant_param_matrix <- matrix(constant_param_vec, nrow = nrow(rsd_vec), 
#                                 ncol=length(constant_param_vec), byrow = TRUE)
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



# constant_env <- outs[[1]]
# variable_env <- outs[[2]]
# 
# # Look at the time-series
# par(mfrow=c(2,1), mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1), las=1)
# matplot(constant_env[,1:4], type="l", xlab="Growing Season", ylab="Biomass", main="Constant Environment", lty=1)
# matplot(variable_env[,1:4], type="l", xlab="Growing Season", ylab="Biomass", main="Fluctuating Environment", lty=1)



####
####  PLOTS FOR INFOGRAPHIC
####
library(ggplot2)
library(viridis)
library(ggthemes)

for(i in 1:length(outs)){
  dosim <- as.data.frame(outs[[i]][4501:5000,5:8])
  colnames(dosim) <- c("Species 1","Species 2","Species 3","Species 4")
  dosim$total <- rowSums(dosim)
  colnames(dosim)[5] <- "Total Biomass"
  dosim$iteration <- 1:500
  dosim_long <- melt(dosim, id.vars = "iteration")
  cv <- round(sd(dosim$`Total Biomass`) / mean(dosim$`Total Biomass`),2)
  # title <- bquote("CV =" ~ .(cv) ~ "and" ~ sigma[E] ~ "=" ~ .(sig_e_vec[i]))
  title <- paste("CV =", cv)
  
  ggplot(dosim_long, aes(x=iteration, y=value, color=variable))+
    geom_line()+
    scale_color_viridis(discrete=TRUE, direction = -1, end = 0.9, name="")+
    xlab("Time")+
    ylab("Biomass")+
    theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
          panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
          panel.background = element_rect("white"),
          axis.line.x = element_line("white"),
          axis.line.y = element_line("white"),
          axis.title = element_text(size=16))+
    ggtitle(title)+
    guides(color=FALSE)
  ggsave(paste0("/Users/atredenn/Desktop/var_",i,"_ts.png"), width = 4, height = 2, units = "in", dpi = 72)
  
}

