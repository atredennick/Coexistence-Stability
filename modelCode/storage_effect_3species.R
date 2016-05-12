##  Semi-discrete storage effect model of two species 
##  coexisting on one essential resource.                

##  Authors:  Andrew Tredennick, Peter Adler, and Fred Adler
##  Email:    atredenn@gmail.com
##  Date:     4.24.2015 -- Copied over from relative_nonlinearity.R and 
##                         variable germination included for storage effect.
##            5.21.2016 -- Updated model to have most dynamics occur at the
##                         end of season, discrete time intervals


############################  MODEL DESCRIPTION  ##############################
# Species "split" themselves between a dormant, low-mortaility stage (D) and  #
# a higher mortality, high growth stage (N).                                  #
# The single resource is R.                                                   #
# There are two sources of variability: an environmental cue that drives the  #
# storage effect, and resource variability.                                   #
###############################################################################

# Clear the Workspace
rm(list=ls())


####
#### Initial Conditions, Global Variables, and Parameters ----------------------
####
seasons <- 1000                   # number of seasons to simulate
seasons_to_exclude <- 100        # initial seasons to exclude from plots
days_to_track <- 20              # number of days to simulate in odSolve
DNR <- c(D=c(1,1,1),N=c(1,1,1),R=20) # initial conditions
Rmu <- 3                         # mean resource pulse (on log scale)
Rsd_annual <- 0                  # std dev of resource pulses (on log scale)
sigE <- 0.01                        # environmental cue variance
rho12 <- -1                         # environmental cue correlation between species
rho13 <- -1
rho23 <- -1

# Within-season parameters
parms <- list(
  r = c(5,5,5),                    # max growth rate for each species
  a = c(5,5,5),                    # rate parameter for Hill function 
  b = c(20,20,20),                  # shape parameter for Hill function
  eps = c(0.5,0.5,0.5)               # resource-to-biomass efficiency
)

# parms <- list (
#   r = c(1,5),                    # max growth rate for each species
#   a = c(2,5),                    # rate parameter for Hill function 
#   b = c(2.5,20),                 # shape parameter for Hill function
#   eps = c(0.5,0.5)               # resource-to-biomass efficiency
# )

# End-of-season transition parameters
alpha1 <- 0.50                   # live-to-dormant biomass fraction; spp1
alpha2 <- 0.49                   # live-to-dormant biomass fraction; spp2
alpha3 <- 0.48
beta1 <- 0                       # adult survivorship; spp1 (0 if annual, >0 if perennial)
beta2 <- 0                       # adult survivorship; spp2 (0 if annual, >0 if perennial)
beta3 <- 0
eta1 <- 0.1                      # dormant mortality; spp1
eta2 <- 0.1                      # dormant mortality; spp2
eta3 <- 0.1
theta1 <- 0                      # resource recycling fraction; spp1
theta2 <- 0                      # resource recycling fraction; spp2
theta3 <- 0
nu <- 0                          # resource carry-over fraction

####
#### Load relevant libraries ----------------------------------------
####
library('deSolve')
library('mvtnorm')



####
#### Model functions -------------------------------------------------
####
## Continuous model
updateNR <- function(t, NR, parms){
  with(as.list(c(NR, parms)), {
    dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
    dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
    dN3dt = N3*eps[3]*(uptake_R(r[3], R, a[3], b[3]))
    dRdt = -1 * (dN1dt/eps[1] + dN2dt/eps[2] + dN3dt/eps[3])
    list(c(dN1dt, dN2dt, dN3dt, dRdt)) #output
  })
}

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, a, b){
  return((r*R^a) / (b^a + R^a))
}

## Discrete model
update_DNR <- function(t,DNR,gammas,alpha1,alpha2,alpha3,eta1,eta2,eta3,beta1,beta2,beta3,theta1,theta2,theta3,nu){
  with (as.list(DNR),{
    g1 <- gammas[1]
    g2 <- gammas[2]
    g3 <- gammas[3]
    D1new <- alpha1*N1 + D1*(1-g1)*(1-eta1)
    D2new <- alpha2*N2 + D2*(1-g2)*(1-eta2)
    D3new <- alpha3*N3 + D3*(1-g3)*(1-eta3)
    N1new <- beta1*(1-alpha1)*N1 + g1*(D1+(alpha1*N1))*(1-eta1)
    N2new <- beta2*(1-alpha2)*N2 + g2*(D2+(alpha2*N2))*(1-eta2)
    N3new <- beta3*(1-alpha3)*N3 + g3*(D3+(alpha3*N3))*(1-eta3)
    Rnew <- theta1*(1-alpha1)*N1 + theta2*(1-alpha2)*N2 + theta3*(1-alpha3)*N3 + nu*R + Rvector[t]
    return(c(D1new, D2new, D3new, N1new, N2new, N3new, Rnew))
  })
}


####
#### Simulate model -----------------------------------------------------
####
days <- c(1:days_to_track)

# Get "germination" fractions for each year
getG <- function(sigE, rho12, rho13, rho23, nTime){
  varcov <- matrix(c(sigE, rho12*sigE, rho13*sigE, rho12*sigE, sigE, rho23*sigE, rho13*sigE, rho23*sigE, sigE), 3, 3, byrow = TRUE)
  varcov <- Matrix::nearPD(varcov)$mat
  e <- MASS::mvrnorm(n = nTime, mu = c(0,0,0), Sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}

# Loop over seasons
nmsDNR <- names(DNR)
NR <- DNR[c("N1", "N2", "N3", "R")] 
nmsNR <- names(NR)
gVec <- getG(sigE = sigE, rho12 = rho12, rho13=rho13, rho23=rho23, nTime = seasons)
Rvector <- rlnorm(seasons, Rmu, Rsd_annual)
saved_outs <- matrix(ncol=7, nrow=seasons)
ptm <- proc.time() # start time
for(season_now in 1:seasons){
  output <- ode(y = NR, times=days,
                func = updateNR, parms = parms)
  NR <- output[nrow(output),nmsNR]
  DNR <- c(DNR[c("D1","D2", "D3")], NR)
  saved_outs[season_now,] <- DNR
  names(DNR) <- nmsDNR
  DNR <- update_DNR(season_now, DNR, gVec[season_now,],
                    alpha1=alpha1,alpha2=alpha2, alpha3=alpha3,
                    eta1=eta1,eta2=eta2, eta3=eta3,
                    beta1=beta1,beta2=beta2, beta3=beta3,
                    theta1=theta1,theta2=theta2, theta3=theta3,nu=nu)
  names(DNR) <- nmsDNR
  NR <- DNR[c("N1", "N2", "N3", "R")] 
  names(NR) <- nmsNR
}
proc.time() - ptm

cols2plot <- c(4,5,6)
par(mfrow=c(1,1),mar=c(4,4,1,1),tcl=-0.2,mgp=c(2,0.5,0))
matplot(saved_outs[(seasons_to_exclude+1):seasons,cols2plot], type="l", 
        xlab="Season", ylab="Abundance")
# lines(rowSums(saved_outs[101:seasons,cols2plot]), type="l", col="blue", lwd=2)

tot_biomass <- rowSums(saved_outs[,cols2plot])
sd(tot_biomass)
mean(tot_biomass)
sd(tot_biomass)/mean(tot_biomass)
 