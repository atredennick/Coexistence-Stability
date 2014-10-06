####
#### Storage Effect model (annual plants)
####

# Following Adler and Drake, AmNat (2008): http://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1203&context=wild_facpub

library(mvtnorm)

####
#### PARAMETERS --------------------------
####

sigE <- c(0,0.4,1,2.5,5,7.5,10)
rho <- c(0.5,0,-0.5)
s <- c(0.5, 0.5)
alpha <- c(1,1)
lambda <- c(101,99)

nTime <- 500


####
#### Germination function -----------------
####

getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}


####
#### Update population function -----------
####

updateN <- function(g, s, alpha, lambda, lastN){
  newN <- numeric(2)
  for(i in 1:2){
    newN[i] <- lastN[i]*s[i]*(1-g[i]) + ((lambda[i]*g[i]*lastN[i]) / (1 + (alpha[i]*g[i]*lastN[i] + alpha[-i]*g[-i]*lastN[-i])))
    newN[i] <- rpois(1, newN[i])
  }
  return(newN)
}


####
#### Simulations ---------------------------
####

gSeries <- getG(sigE[4], rho[1], nTime)
N <- matrix(data = NA, nrow = nTime, ncol = 2)
N[1,] <- c(100,100)
for(t in 2:nTime){
  N[t,] <- updateN(g=gSeries[t,], s, alpha, lambda, lastN=N[t-1,])
}

matplot(c(1:nTime), N, type="l", lwd=2)



