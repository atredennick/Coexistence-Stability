##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##            

##  These are the core model functions for the multi-species case
##  that can be called from simulation wrappers.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Last update: 6-06-2016


simulate_model <- function(seasons, days_to_track, Rmu, 
                           Rsd_annual, sigE, rho, 
                           alphas, betas, etas, thetas, rs,
                           as, bs, epss, nu,
                           Ds, Ns, R, varcov) {
  require('deSolve') # for solving continuous differential equations
  require('mvtnorm') # for multivariate normal distribution functions
  
  nspp <- length(alphas)
  
  ##  Assign parameter values to appropriate lists
  DNR = c(D=Ds, N=Ns, R=R) # initial conditions
  parms <- list (
    r = rs,       # max growth rate for each species
    a = as,       # rate parameter for Hill function 
    b = bs,       # shape parameter for Hill function
    eps = epss    # resource-to-biomass efficiency
  )
  
  
  ####
  ####  Sub-Model functions
  ####
  
  ########################
  ######  2 SPECIES ######
  ########################
  if(nspp == 2) {
    
    ## Continuous model
    updateNR <- function(t, NR, parms) {
      with(as.list(c(NR, parms)), {
        dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
        dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
        dRdt = -1 * (dN1dt/eps[1] + dN2dt/eps[2])
        list(c(dN1dt, dN2dt, dRdt)) #output
      })
    } # end continuous function
    
    ## Discrete model
    update_DNR <- function(t,DNR,gammas,alphas,etas,
                           betas,thetas,nu) {
      with (as.list(DNR),{
        g1 <- gammas[1]
        g2 <- gammas[2]
        D1new <- alphas[1]*N1 + D1*(1-g1)*(1-etas[1])
        D2new <- alphas[2]*N2 + D2*(1-g2)*(1-etas[2])
        N1new <- betas[1]*(1-alphas[1])*N1 + g1*(D1+(alphas[1]*N1))*(1-etas[1])
        N2new <- betas[2]*(1-alphas[2])*N2 + g2*(D2+(alphas[2]*N2))*(1-etas[2])
        Rnew <- thetas[1]*(1-alphas[1])*N1 + thetas[2]*(1-alphas[2])*N2 + nu*R + Rvector[t]
        return(c(D1new, D2new, N1new, N2new, Rnew))
      })
    } # end discrete function
    
  }
  ############################
  ######  END 2 SPECIES ######
  ############################
  
  
  
  ########################
  ######  3 SPECIES ######
  ########################
  if(nspp == 3) {
    
    ## Continuous model
    updateNR <- function(t, NR, parms){
      with(as.list(c(NR, parms)), {
        dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
        dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
        dN3dt = N3*eps[3]*(uptake_R(r[3], R, a[3], b[3]))
        dRdt = -1 * (dN1dt/eps[1] + dN2dt/eps[2] + dN3dt/eps[3])
        list(c(dN1dt, dN2dt, dN3dt, dRdt)) #output
      })
    } # end continuous function
    
    ## Discrete model
    update_DNR <- function(t,DNR,gammas,alphas,etas,
                           betas,thetas,nu) {
      with (as.list(DNR),{
        g1 <- gammas[1]
        g2 <- gammas[2]
        g3 <- gammas[3]
        D1new <- alphas[1]*N1 + D1*(1-g1)*(1-etas[1])
        D2new <- alphas[2]*N2 + D2*(1-g2)*(1-etas[2])
        D3new <- alphas[3]*N3 + D3*(1-g3)*(1-etas[3])
        N1new <- betas[1]*(1-alphas[1])*N1 + g1*(D1+(alphas[1]*N1))*(1-etas[1])
        N2new <- betas[2]*(1-alphas[2])*N2 + g2*(D2+(alphas[2]*N2))*(1-etas[2])
        N3new <- betas[3]*(1-alphas[3])*N3 + g3*(D3+(alphas[3]*N3))*(1-etas[3])
        Rnew <- thetas[1]*(1-alphas[1])*N1 + thetas[2]*(1-alphas[2])*N2 + thetas[3]*(1-alphas[3])*N3 + nu*R + Rvector[t]
        return(c(D1new, D2new, D3new, N1new, N2new, N3new, Rnew))
      })
    } # end discrete function
    
  }
  ############################
  ######  END 3 SPECIES ######
  ############################
  
  
  
  ########################
  ######  4 SPECIES ######
  ########################
  if(nspp == 4) {
    
    ## Continuous model
    updateNR <- function(t, NR, parms){
      with(as.list(c(NR, parms)), {
        dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
        dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
        dN3dt = N3*eps[3]*(uptake_R(r[3], R, a[3], b[3]))
        dN4dt = N4*eps[4]*(uptake_R(r[4], R, a[4], b[4]))
        dRdt = -1 * (dN1dt/eps[1] + dN2dt/eps[2] + dN3dt/eps[3] + dN4dt/eps[4])
        list(c(dN1dt, dN2dt, dN3dt, dN4dt, dRdt)) #output
      })
    } # end continuous function
    
    ## Discrete model
    update_DNR <- function(t,DNR,gammas,alphas,etas,
                           betas,thetas,nu) {
      with (as.list(DNR),{
        g1 <- gammas[1]
        g2 <- gammas[2]
        g3 <- gammas[3]
        g4 <- gammas[4]
        D1new <- alphas[1]*N1 + D1*(1-g1)*(1-etas[1])
        D2new <- alphas[2]*N2 + D2*(1-g2)*(1-etas[2])
        D3new <- alphas[3]*N3 + D3*(1-g3)*(1-etas[3])
        D4new <- alphas[4]*N4 + D4*(1-g4)*(1-etas[4])
        N1new <- betas[1]*(1-alphas[1])*N1 + g1*(D1+(alphas[1]*N1))*(1-etas[1])
        N2new <- betas[2]*(1-alphas[2])*N2 + g2*(D2+(alphas[2]*N2))*(1-etas[2])
        N3new <- betas[3]*(1-alphas[3])*N3 + g3*(D3+(alphas[3]*N3))*(1-etas[3])
        N4new <- betas[4]*(1-alphas[4])*N4 + g4*(D4+(alphas[4]*N4))*(1-etas[4])
        Rnew <- thetas[1]*(1-alphas[1])*N1 + thetas[2]*(1-alphas[2])*N2 + thetas[3]*(1-alphas[3])*N3 + thetas[4]*(1-alphas[4])*N4 + nu*R + Rvector[t]
        return(c(D1new, D2new, D3new, D4new, N1new, N2new, N3new, N4new, Rnew))
      })
    } # end discrete function
    
  }
  ############################
  ######  END 4 SPECIES ######
  ############################
  
  
  
  ##  Resource uptake function (Hill function)
  uptake_R <- function(r, R, a, b) {
    return((r*R^a) / (b^a + R^a))
  }
  
  
  
  
  ####
  #### Simulate model -----------------------------------------------------
  ####
  days <- c(1:days_to_track)
  
  # Get "germination" fractions for each year
  getG <- function(sigE, varcov, nTime){
    e <- MASS::mvrnorm(n = nTime, mu = rep(0,nspp), Sigma = varcov)
    g <- exp(e) / (1+exp(e))
    return(g)
  }
  
  # Loop over seasons
  nmsDNR <- names(DNR)
  dormants <- grep("D", names(DNR))
  NR <- DNR[-dormants] 
  nmsNR <- names(NR)
  gVec <- getG(sigE = sigE, varcov, nTime = seasons)
  Rvector <- rlnorm(seasons, Rmu, Rsd_annual)
  saved_outs <- matrix(ncol=length(DNR), nrow=seasons+1)
  saved_outs[1,] <- DNR 

  for(season_now in 1:seasons) {
    output <- ode(y = NR, times=days,
                  func = updateNR, parms = parms)
    NR <- output[nrow(output),nmsNR]
    dormants <- grep("D", names(DNR))
    DNR <- c(DNR[dormants], NR)
    saved_outs[season_now+1,] <- DNR
    names(DNR) <- nmsDNR
    DNR <- update_DNR(season_now, DNR, gVec[season_now,],
                      alphas=alphas, etas=etas,
                      betas=betas, thetas=thetas, nu=nu)
    names(DNR) <- nmsDNR
    NR <- DNR[-dormants] 
    names(NR) <- nmsNR
  } # next season
  
  return(saved_outs)
  
} #end simulation function

