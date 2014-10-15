# coexistence via relative nonlinearity
# two species competing for a shared, fluctuating resource
rm(list=ls())

####
#### Parameters
####

# initial conditions
maxDays <- 40 #Rfreq*cycles
years <- 500
burn.in <- 120
NR <- c(N=c(1,1),R=50)
Rmean <- 50

parmsB <- list(
  a = c(0.23,0.2),  # max growth rate for spp 1 and 2
  b1 = c(25,10),  # right offset for growth rates 
  b2 = c(0.08,0.45)  # rates at which max is approached
)


#---------------------------
# define functions
#---------------------------
library(deSolve)
updateNR = function(t, NR, parms) { 
  #### RRmod "Resource Ratio Model" 
  #### Continuous two species resource model with an essential resource
  with(as.list(c(NR, parms)), { 
    dN1dt = N1*b[1]*(c[1]*w[1]*R - m[1]) ##spp 1 growth
    dN2dt = N2*b[2]*(c[2]*w[2]*R - m[2]) ##spp 2 growth
    dRdt = -1 * (c[1]*N1 + c[2]*N2) ##resource draw down
    list(c(dN1dt, dN2dt, dRdt))
  }) 
}

getB <- function(parmsB, R){
  #   attach(parms)
  spp1 <- parmsB$a[1]*exp(-parmsB$b1[1]*exp(-parmsB$b2[1]*R))
  spp2 <- parmsB$a[2]*exp(-parmsB$b1[2]*exp(-parmsB$b2[2]*R))   ## sp. 2 biomass assimilation
  #   detach(parms)
  return(c(spp1, spp2))
}

####
#### Run simulations
####
Nstart <- 5
Nsave <- matrix(nrow = years, ncol=2)
Nsave[1,] <- c(Nstart,Nstart)

#get vector of means for resource
# set.seed(1653)
muVec <- rnorm(years, 0, .2)

for(g in 2:years){
  days = 1:maxDays
  mu <- muVec[g]
  NR[3] <- rnorm(1, Rmean, 13)
  # species parameters 
  cNow <- getB(parmsB, NR[3])
  parms <- list(
    b = c(0.1, 0.1), #biomass conversion efficiency
    c = as.numeric(cNow), #consumption rates of resource
    w = c(0.5, 0.5), #biomass assimilation rate
    m = c(5, 5) #loss of assimilated biomass
  )
  output = as.data.frame(ode(y = NR, times = days, func = updateNR, parms = parms))
  N = output[,2:3]
  nMax <- c(max(N[1]), max(N[2]))
  seedFrac1 <- nMax[1] / (sum(nMax))
  NR <- c(N=c(seedFrac1*(Nstart*2), (1-seedFrac1)*(Nstart*2)),R=NA)
  Nsave[g,] <- nMax
}

par(mfrow=c(1,3),mar=c(4,4,1,1))
matplot(c(1:years),Nsave,type="l",xlab="Time",ylab="N")

gomp <- function(parmsB,R){
  #   attach(parms)
  out1 <- parmsB$a[1]*exp(-parmsB$b1[1]*exp(-parmsB$b2[1]*R))
  out2 <- parmsB$a[2]*exp(-parmsB$b1[2]*exp(-parmsB$b2[2]*R))
  #   detach(parms)
  return(cbind(out1,out2))
}
R <- 0:120
tmp<-gomp(parmsB=parmsB,R)
matplot(R,tmp,type="l",ylab="Resource uptake")

matplot(output[,1], output[,2:3], type="l")

