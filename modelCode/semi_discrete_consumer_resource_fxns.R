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
## Continuous model
updateDNR <- function(t, DNR, parms){
  with(as.list(c(DNR, parms)), {
    dD1dt = N1*a[1] - (mD[1]*D1)
    dD2dt = N2*a[2] - (mD[2]*D2)
    dN1dt = N1*(uptake_R(r[1], R, alpha[1], beta[1]) - mN[1] - a[1])
    dN2dt = N2*(uptake_R(r[2], R, alpha[2], beta[2]) - mN[2] - a[2])
    dRdt = -1 * ((dN1dt + mN[1]*N1) + (dN2dt + mN[2]*N2))
    list(c(dD1dt, dD2dt, dN1dt, dN2dt, dRdt)) #output
  })
}

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, alpha, beta){
  return((r*R^alpha) / (beta^alpha + R^alpha))
}

## Discrete model
update_DNR <- function(t,DNR,gs){
  with (as.list(DNR),{
    g1 <- gs[1]
    g2 <- gs[2]
    N1 <- D1*g1
    N2 <- D2*g2
    D1 <- D1*(1-g1)
    D2 <- D2*(1-g2)
    R <- Rvector[t]
    return(c(D1, D2, N1, N2, R))
  })
}


# Get "germination" fractions for each year
getG <- function(sigE, rho, nTime){
  varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
  e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}

