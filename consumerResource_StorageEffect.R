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
b <- c(0.5, 0.5)
c <- c(0.1, 0.1)
w <- c(0.2, 0.2)
m <- c(0.15, 0.15)
d <- 0.5
nTime <- 1000


####
#### New biomass assimilation function ------------------------
####

getB <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  b <- exp(e) / (1+exp(e))
  return(b)
}


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
  newR <- max(0, s - d*lastR - sum(consumption))
  return(list(newN, newR))
}


####
#### Model simulation
####

bSeries <- getB(sigE[4], rho[3], nTime)
N <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,] <- c(1,1)
R <- numeric(nTime)
R[1] <- 0.1

for(t in 2:nTime){
  s <- rtruncnorm(n = 1, a = 0, b = Inf, mean = 15, sd = 5)
  thisStep <- updateN(b=bSeries[t,], c, w, s, d, lastN=N[t-1,], lastR=R[t-1])
#   N[t,] <- rpois(2, thisStep[[1]])
  N[t,] <- thisStep[[1]]
  R[t] <- thisStep[[2]]
}

par(mfrow=c(3,1))
matplot(c(1:nTime), N, type="l", lwd=2)
plot(c(1:nTime), R, type="l", lwd=1)
matplot(c(1:nTime), bSeries, type="l", lwd=2)



