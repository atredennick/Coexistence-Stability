##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##            

##  These are the core model functions that can be called 
##    from simulation wrappers.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Last update: 4-27-2016


simulate_model <- function(seasons, days_to_track, Rmu, 
                           Rsd_annual, sigE, rho, 
                           alpha1, alpha2, beta1, beta2,
                           eta1,eta2, theta1, theta2, nu, r1,
                           r2, a1, a2, b1, b2, eps1, eps2,
                           D1, D2, N1, N2, R) {
  require('deSolve') # for solving continuous differential equations
  require('mvtnorm') # for multivariate normal distribution functions
  
  ##  Assign parameter values to appropriate lists
  DNR = c(D=c(D1,D2),N=c(N1,N2),R=R) # initial conditions
  parms <- list (
    r = c(r1,r2),         # max growth rate for each species
    a = c(a1,a2),         # rate parameter for Hill function 
    b = c(b1,b2),         # shape parameter for Hill function
    eps = c(eps1,eps2)    # resource-to-biomass efficiency
  )
  
  
  ####
  ####  Sub-Model functions
  ####
  ## Continuous model
  updateNR <- function(t, NR, parms){
    with(as.list(c(NR, parms)), {
      dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
      dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
      dRdt = -1 * (dN1dt/eps[1] + dN2dt/eps[2])
      list(c(dN1dt, dN2dt, dRdt)) #output
    })
  }
  
  ##  Resource uptake function (Hill function)
  uptake_R <- function(r, R, a, b){
    return((r*R^a) / (b^a + R^a))
  }
  
  ## Discrete model
  update_DNR <- function(t,DNR,gammas,alpha1,alpha2,eta1,eta2,beta1,beta2,theta1,theta2,nu){
    with (as.list(DNR),{
      g1 <- gammas[1]
      g2 <- gammas[2]
      D1 <- alpha1*N1 + D1*(1-g1)*(1-eta1)
      D2 <- alpha2*N2 + D2*(1-g2)*(1-eta2)
      N1 <- beta1*(1-alpha1)*N1 + g1*(D1+(alpha1*N1))*(1-eta1)
      N2 <- beta2*(1-alpha2)*N2 + g2*(D2+(alpha2*N2))*(1-eta2)
      R <- theta1*(1-alpha1)*N1 + theta2*(1-alpha2)*N2 + nu*R + Rvector[t]
      return(c(D1, D2, N1, N2, R))
    })
  }
  
  
  ####
  #### Simulate model -----------------------------------------------------
  ####
  days <- c(1:days_to_track)
  
  # Get "germination" fractions for each year
  getG <- function(sigE, rho, nTime){
    varcov <- matrix(c(sigE, rho*sigE, rho*sigE, sigE), 2, 2)
    e <- rmvnorm(n = nTime, mean = c(0,0), sigma = varcov)
    g <- exp(e) / (1+exp(e))
    return(g)
  }
  
  # Loop over seasons
  nmsDNR <- names(DNR)
  NR <- DNR[c("N1", "N2", "R")] 
  nmsNR <- names(NR)
  gVec <- getG(sigE = sigE, rho = rho, nTime = seasons)
  Rvector <- rlnorm(seasons, Rmu, Rsd_annual)
  saved_outs <- matrix(ncol=5, nrow=seasons)

  for(season_now in 1:seasons) {
    output <- ode(y = NR, times=days,
                  func = updateNR, parms = parms)
    NR <- output[nrow(output),nmsNR]
    DNR <- c(DNR[c("D1","D2")], NR)
    saved_outs[season_now,] <- DNR
    names(DNR) <- nmsDNR
    DNR <- update_DNR(season_now, DNR, gVec[season_now,],
                      alpha1=alpha1,alpha2=alpha2,
                      eta1=eta1,eta2=eta2,
                      beta1=beta1,beta2=beta2,
                      theta1=theta1,theta2=theta2,nu=nu)
    names(DNR) <- nmsDNR
    NR <- DNR[c("N1", "N2", "R")] 
    names(NR) <- nmsNR
  } # next season
  
  return(saved_outs)
  
} #end simulation

