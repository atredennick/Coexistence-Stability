##  Wrapper for running simulations of the storage effect
##    model with variable resource.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

# clear the workspace
rm(list=ls())


####
#### Initial conditions and global variables ------------------------
####
maxTime <- 2500 
burn.in <- 1000
DNR <- c(D=c(1,1),N=c(1,1),R=100)
Rmu <- 2                                #mean resource pulse (on log scale)
Rsd <- seq(0,3,by=0.25)                 #std dev of resource pulses (on log scale)
sigE <- c(0, 0.4, 1, 2.5, 5, 7.5, 10)   #environmental cue variability
rho_vec <- seq(-1,1,by=0.25)            #environmental cue correlation between species
parms <- list(
  r = c(5,5),                           #max growth rate for genotype A and a
  k1 = c(20,20),                        #right offset for growth rates 
  k2 = c(0.08,0.08),                    #rates at which max is approached
  mN = c(0.5,0.5),                      #live biomass loss (mortality) rates 
  mD = c(0.001, 0.001)                  #dormant biomass loss (mortality) rates
)

resource_sims <- length(Rsd)
env_cue_sims <- length(sigE)
rho_sims <- length(rho_vec)
sims_per_level <- 2

####
#### Load relevant libraries and source the model --------------------
####
library('deSolve')
library('mvtnorm')
source("semi_discrete_consumer_resource_fxns.R")


####
####  Loop through simulation sets -----------------------------------
####
out_sims <- data.frame(sim=NA, rho=NA, sig_e=NA, sig_r=NA,
                       sd_n=NA, mu_n=NA, cv_n=NA, buffer=NA)
simTime <- seq(1,maxTime,by=1) #creates a vector of time steps to get output from deSolve
for(resource in 1:resource_sims){
  for(cue in 1:env_cue_sims){
    for(rho in 1:rho_sims){
      for(sim in 1:sims_per_level){
        gVec <- getG(sigE = sigE[cue], rho = rho_vec[rho], nTime = maxTime)
        gVec1 <- gVec[,1]
        gVec2 <- gVec[,2]
        
        # Set random resource fluctuations
        Rvector <- rlnorm(maxTime,Rmu,Rsd[resource])
        
        # Run the model
        output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR, parms = parms,
                                   events = list(func = gfun, times=simTime)))
        
        # Collect outputs
        biomass_sd <- sd(rowSums(output[burn.in:maxTime,2:3])) 
        biomass_mu <- mean(rowSums(output[burn.in:maxTime,2:3]))
        biomass_cv <- biomass_sd/biomass_mu 
        resource_cv <- sd(output[burn.in:maxTime,4])/mean(output[burn.in:maxTime,4])
        buffer <- biomass_cv/resource_cv
        out_sims$sim <- sim
        out_sims$rho <- rho_vec[rho]
        out_sims$sig_e <- sigE[cue]
        out_sims$sig_r <- Rsd[resource]
        out_sims$sd_n <- biomass_sd
        out_sims$mu_n <- biomass_mu
        out_sims$cv_n <- biomass_cv
        out_sims$buffer <- buffer
          
        print(paste("Done with resource", resource, "for sigE", cue, "and rho", rho))
      } #next simulation
    } #next rho
  } #next cue variability level
} #next resource variability level


