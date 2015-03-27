##  Semi-discrete storage effect model of two species   
##    coexisting on one essential resource               

##  These are the core model functions that can be called 
##    from simulation wrappers.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Date:   3.27.2015

####
####  Model functions ------------------------------------
####
##  Continuous model
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = -(mD[1]*D1)
    dD2dt = -(mD[2]*D2)
    dN1dt = N1*(r[1]*exp(-k1[1]*exp(-k2[1]*R)) - mN[1])
    dN2dt = N2*(r[2]*exp(-k1[2]*exp(-k2[2]*R)) - mN[2])
    dRdt = -1* ((dN1dt + mN[1]*N1) + (dN2dt + mN[2]*N2))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt)) #output
  })
}

##  Discrete model
gfun <- function(t, y, parms){
  with (as.list(y),{
    g1 <- gVec1[t]
    g2 <- gVec2[t]
    D1 <- D1 - (N1+D1)*g1 + N1
    D2 <- D2 - (N2+D2)*g2 + N2
    N1 <- 0+(N1+D1)*g1
    N2 <- 0+(N2+D2)*g2
    R <- R + Rvector[t]
    return(c(D1, D2, N1, N2, R))
  })
}

##  Get "germination" fractions for each year
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}

