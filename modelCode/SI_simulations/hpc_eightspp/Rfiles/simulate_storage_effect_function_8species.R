##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "The relationship between species richness and ecosystem variability is    ##
##  shaped by the mechanism of coexistence"                                   ##
##                                                                       2017 ##
##============================================================================##            

##  These are the model functions for the storage effect with many species that
##  can be called from simulation wrappers. This is an 8 species version.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Last update: 2-22-2017


simulate_model <- function(seasons, days_to_track, Rmu, 
                           Rsd_annual, sigE, rho, 
                           alpha1, alpha2, alpha3, alpha4,
                           alpha5, alpha6, alpha7, alpha8,
                           eta1, eta2, eta3, eta4,
                           eta5, eta6, eta7, eta8,
                           r1, r2, r3, r4,
                           r5, r6, r7, r8,
                           a1, a2, a3, a4,
                           a5, a6, a7, a8,
                           b1, b2, b3, b4,
                           b5, b6, b7, b8,
                           eps1, eps2, eps3, eps4,
                           eps5, eps6, eps7, eps8,
                           D1, D2, D3, D4,
                           D5, D6, D7, D8,
                           N1, N2, N3, N4,
                           N5, N6, N7, N8, R) {
  
  require('deSolve') # for solving continuous differential equations
  require('mvtnorm') # for multivariate normal distribution functions
  
  ##  Assign parameter values to appropriate lists
  DNR <- c(D=c(D1,D2,D3,D4,D5,D6,D7,D8),   # initial dormant state abundance
           N=c(N1,N2,N3,N4,N5,N6,N7,N8),   # initial live state abundance
           R=R)                            # initial resource level
  
  parms <- list (
    r   = c(r1,r2,r3,r4,r5,r6,r7,r8),   # max growth rate for each species
    a   = c(a1,a2,a3,a4,a5,a6,a7,a8),   # rate parameter for Hill function 
    b   = c(b1,b2,b3,b4,b5,b6,b7,b8),   # shape parameter for Hill function
    eps = c(eps1,eps2,eps3,eps4,eps5,eps6,eps7,eps8) # resource-to-biomass efficiency
  )
  
  
  ####
  ####  Sub-Model functions ----------------------------------------------------
  ####
  ## Continuous model
  updateNR <- function(t, NR, parms){
    with(as.list(c(NR, parms)), {
      dN1dt = N1*eps[1]*(uptake_R(r[1], R, a[1], b[1]))
      dN2dt = N2*eps[2]*(uptake_R(r[2], R, a[2], b[2]))
      dN3dt = N3*eps[3]*(uptake_R(r[3], R, a[3], b[3]))
      dN4dt = N4*eps[4]*(uptake_R(r[4], R, a[4], b[4]))
      dN5dt = N5*eps[5]*(uptake_R(r[5], R, a[5], b[5]))
      dN6dt = N6*eps[6]*(uptake_R(r[6], R, a[6], b[6]))
      dN7dt = N7*eps[7]*(uptake_R(r[7], R, a[7], b[7]))
      dN8dt = N8*eps[8]*(uptake_R(r[8], R, a[8], b[8]))
      dRdt  = -1 * (dN1dt/eps[1] + dN2dt/eps[2] + dN3dt/eps[3] + dN4dt/eps[4] +
                    dN5dt/eps[5] + dN6dt/eps[6] + dN7dt/eps[7] + dN8dt/eps[8])
      list(c(dN1dt, dN2dt, dN3dt, dN4dt, dN5dt, dN6dt, dN7dt, dN8dt, dRdt)) # output as list
    })
  } # end continuous function
  
  ## Discrete model
  update_DNR <- function(t, DNR, gammas,
                         alpha1, alpha2, alpha3, alpha4,
                         alpha5, alpha6, alpha7, alpha8,
                         eta1, eta2, eta3, eta4,
                         eta5, eta6, eta7, eta8) {
    with (as.list(DNR),{
      g1    <- gammas[1]
      g2    <- gammas[2]
      g3    <- gammas[3]
      g4    <- gammas[4]
      g5    <- gammas[5]
      g6    <- gammas[6]
      g7    <- gammas[7]
      g8    <- gammas[8]
      D1new <- (1-g1)*(alpha1*N1 + D1)*(1-eta1)
      D2new <- (1-g2)*(alpha2*N2 + D2)*(1-eta2)
      D3new <- (1-g3)*(alpha3*N3 + D3)*(1-eta3)
      D4new <- (1-g4)*(alpha4*N4 + D4)*(1-eta4)
      D5new <- (1-g5)*(alpha5*N5 + D5)*(1-eta5)
      D6new <- (1-g6)*(alpha6*N6 + D6)*(1-eta6)
      D7new <- (1-g7)*(alpha7*N7 + D7)*(1-eta7)
      D8new <- (1-g8)*(alpha8*N8 + D8)*(1-eta8)
      N1new <- g1*(alpha1*N1 + D1)*(1-eta1)
      N2new <- g2*(alpha2*N2 + D2)*(1-eta2)
      N3new <- g3*(alpha3*N3 + D3)*(1-eta3)
      N4new <- g4*(alpha4*N4 + D4)*(1-eta4)
      N5new <- g5*(alpha5*N5 + D5)*(1-eta5)
      N6new <- g6*(alpha6*N6 + D6)*(1-eta6)
      N7new <- g7*(alpha7*N7 + D7)*(1-eta7)
      N8new <- g8*(alpha8*N8 + D8)*(1-eta8)
      Rnew  <- Rvector[t]
      return(c(D1new, D2new, D3new, D4new, D5new, D6new, D7new, D8new,
               N1new, N2new, N3new, N4new, N6new, N6new, N7new, N8new, Rnew))
    })
  }
  
  ##  Resource uptake function (Hill function)
  uptake_R <- function(r, R, a, b) {
    return((r*R^a) / (b^a + R^a))
  }
  
  ##  Generate germination fractions
  getG <- function(sigE, rho, nTime, num_spp) {
    varcov       <- matrix(rep(rho*sigE,num_spp*2), num_spp, num_spp)
    diag(varcov) <- sigE
    # if(sigE > 0) { varcov <- Matrix::nearPD(varcov)$mat } # crank through nearPD to fix rounding errors 
    # varcov <- as.matrix(varcov)
    e      <- rmvnorm(n = nTime, mean = rep(0,num_spp), sigma = varcov)
    g      <- exp(e) / (1+exp(e))
    return(g)
  }
  
  
  ####
  #### Simulate model -----------------------------------------------------
  ####
  days           <- c(1:days_to_track)
  num_spp        <- length(parms$r)
  nmsDNR         <- names(DNR)
  dormants       <- grep("D", names(DNR))
  NR             <- DNR[-dormants] 
  nmsNR          <- names(NR)
  gVec           <- getG(sigE = sigE, rho = rho, nTime = seasons, num_spp = num_spp)
  Rvector        <- rlnorm(seasons, Rmu, Rsd_annual)
  saved_outs     <- matrix(ncol=length(DNR), nrow=seasons+1)
  saved_outs[1,] <- DNR 

  ##  Loop over seasons
  for(season_now in 1:seasons) {
    # Simulate continuous growing  season
    output   <- ode(y = NR, times=days, func = updateNR, parms = parms)
    NR       <- output[nrow(output),nmsNR]
    dormants <- grep("D", names(DNR)) 
    DNR      <- c(DNR[dormants], NR)
    
    # Save end of season biomasses, before discrete transitions
    saved_outs[season_now+1,] <- DNR
    
    names(DNR) <- nmsDNR
    DNR <- update_DNR(season_now, DNR, gVec[season_now,],
                      alpha1 = alpha1, alpha2 = alpha2, 
                      alpha3 = alpha3, alpha4 = alpha4,
                      alpha5 = alpha5, alpha6 = alpha6, 
                      alpha7 = alpha7, alpha8 = alpha8,
                      eta1 = eta1, eta2 = eta2, eta3 = eta3, eta4 = eta4,
                      eta5 = eta5, eta6 = eta6, eta7 = eta7, eta8 = eta8)
    
    names(DNR) <- nmsDNR
    NR         <- DNR[-dormants]  
    names(NR)  <- nmsNR
  } # next season
  
  return(saved_outs)
  
} #end simulation function

