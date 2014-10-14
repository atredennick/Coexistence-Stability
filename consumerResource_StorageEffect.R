####
#### CONSUMER-RESOURCE STORAGE EFFECT MODEL (2 species)
####

#### A. Tredennick // 10-06-2014

library(mvtnorm)
library(truncnorm)

####
#### Stationary model parameters ------------------------------
####

sigE <- c(0,0.4,1,2.5,5,7.5,10)
rho <- c(0.5,0,-0.5)
# b <- c(0.5, 0.5)
c <- c(0.1, 0.1)
w <- c(0.2, 0.2)
m <- c(0.6, 0.6)
d <- 0.5
nTime <- 1000


####
#### New biomass assimilation function ------------------------
####

# getB <- function(sigE, rho, nTime){
#   varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
#   e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
#   b <- exp(e) / (1+exp(e))
#   return(b)
# }
getB <- function(parms, R){
#   attach(parms)
  spp1 <- parms$a[1]*exp(-parms$b1[1]*exp(-parms$b2[1]*R))
  spp2 <- parms$a[2]*exp(-parms$b1[2]*exp(-parms$b2[2]*R))   ## sp. 2 biomass assimilation
#   detach(parms)
  return(c(spp1, spp2))
}

parms <- list(
  a = c(0.4,0.2),  # max growth rate for spp 1 and 2
  b1 = c(1.5,0.8),  # right offset for growth rates 
  b2 = c(0.1,0.5)  # rates at which max is approached
)

gomp <- function(parms,R){
#   attach(parms)
  out1 <- parms$a[1]*exp(-parms$b1[1]*exp(-parms$b2[1]*R))
  out2 <- parms$a[2]*exp(-parms$b1[2]*exp(-parms$b2[2]*R))
#   detach(parms)
  return(cbind(out1,out2))
}
R <- 0:25
tmp<-gomp(parms=parms,R)
par(mfrow=c(1,1),mar=c(4,4,1,1))
matplot(R,tmp,type="l",ylab="Resource uptake")


####
#### Update population and resource function ------------------
####

updateN <- function(b, c, w, s, d, lastN, lastR){
  newN <- numeric(2)
  consumption <- numeric(2)
  for(i in 1:2){
    newN[i] <- lastN[i] + lastN[i]*b[i]*(c[i]*w[i]*lastR - m[i])
    consumption[i] <- c[i]*lastR*lastN[i]
  }
  newR <- max(0, s)
  return(list(newN, newR))
}


####
#### Model simulation
####

# bSeries <- getB(sigE[4], rho[3], nTime)
N <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,] <- c(1,1)
R <- numeric(nTime)
R[1] <- 0.1

for(t in 2:nTime){
#   s <- rtruncnorm(n = 1, a = 0, b = 13*2, mean = 13, sd = 5)
  s <- rlnorm(n=1, mean=4.1, 0.1)
  bTmp <- getB(parms, R[t-1])
  thisStep <- updateN(b=bTmp, c, w, s, d, lastN=N[t-1,], lastR=R[t-1])
#   thisStep <- updateN(b=bSeries[t,], c, w, s, d, lastN=N[t-1,], lastR=R[t-1])
#   N[t,] <- rpois(2, thisStep[[1]])
  N[t,] <- thisStep[[1]]
  R[t] <- thisStep[[2]]
}

par(mfrow=c(2,1))
matplot(c(1:nTime), N, type="l", lwd=2)
plot(c(1:nTime), R, type="l", lwd=1)
# matplot(c(1:nTime), bSeries, type="l", lwd=2)



