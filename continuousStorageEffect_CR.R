#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Continuous-time temporal storage effect model of two species coexisting on one essential resource ##
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# species "split" themselves between a dormant, low-mortaility stage (D) and a higher mortality, high growth stage (N)
# the single resource is R

rm(list=ls())
library(deSolve)

####
#### Parameters
####
maxTime <- 100
c <- c(0.5,0.5)
b <- c(0.2, 0.2)
g <- c(0.2, 0.2)
mD <- c(0.0001, 0.0001)
r <- c(2, 2)
K <- c(1, 1)
mN <- c(0.1, 0.1)
a <- 0.5
S <- 5

####
#### Model function
####
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = c[1]*b[1]*N1 - g[1]*D1 - mD[1]*D1
    dD2dt = c[2]*b[2]*N2 - g[2]*D2 - mD[2]*D2
    dN1dt = N1*(r[1]*R/(K[1]+R)) + g[1]*D1 - c[1]*b[1]*N1 - mN[1]*N1
    dN2dt = N2*(r[2]*R/(K[2]+R)) + g[2]*D2 - c[2]*b[2]*N2 - mN[2]*N2
    dRdt = a*(S-R) - (N1*(r[1]*R/(K[1]+R)) + N2*(r[2]*R/(K[2]+R)))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt))
  })
}

####
#### Simulate model
####
simTime <- c(1:maxTime)
parms <- list(
  c = c,
  b = b,
  g = g,
  mD = mD,
  r = r,
  K = K,
  mN = mN, 
  a = a,
  S = S
)
DNR <- c(D=c(1,1), N=c(1,1),R=50)
output = as.data.frame(ode(y = DNR, times = simTime, func = updateDNR, parms = parms))


####
#### Plot results
####
par(mfrow=c(1,2))
matplot(simTime, output[,4:5], type="l", main="Live Biomass")
matplot(simTime, output[,2:3], type="l", main="Dormant Biomass")


