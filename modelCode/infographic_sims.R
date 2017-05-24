#########################################################################
##  test_runs.R: does quick simulations of the four species model to   ##
##  test parameter combinations                                        ##
#########################################################################

rm(list=ls())                    # Erase the memory
library('mvtnorm')

####
####  THE MODEL FUNCTION -------------------------------------------
####
simulate_model <- function(seasons, days_to_track, Rmu, 
                           Rsd_annual, sigE, rho, 
                           alpha1, alpha2, alpha3, alpha4,
                           eta1, eta2, eta3, eta4,
                           r1, r2, r3, r4,
                           a1, a2, a3, a4,
                           b1, b2, b3, b4,
                           eps1, eps2, eps3, eps4,
                           D1, D2, D3, D4,
                           N1, N2, N3, N4, R, gVec) {
  
  require('deSolve') # for solving continuous differential equations
  require('mvtnorm') # for multivariate normal distribution functions
  
  ##  Assign parameter values to appropriate lists
  DNR <- c(D=c(D1,D2,D3,D4),    # initial dormant state abundance
           N=c(N1,N2,N3,N4),    # initial live state abundance
           R=R)                 # initial resource level
  
  parms <- list (
    r = c(r1,r2,r3,r4),          # max growth rate for each species
    a = c(a1,a2,a3,a4),          # rate parameter for Hill function 
    b = c(b1,b2,b3,b4),          # shape parameter for Hill function
    eps = c(eps1,eps2,eps3,eps4) # resource-to-biomass efficiency
  )
  
  
  ####
  ####  Sub-Model functions
  ####
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
  update_DNR <- function(t, DNR, gammas,
                         alpha1, alpha2, alpha3, alpha4,
                         eta1, eta2, eta3, eta4) {
    with (as.list(DNR),{
      g1    <- gammas[1]
      g2    <- gammas[2]
      g3    <- gammas[3]
      g4    <- gammas[4]
      D1new <- (1-g1)*(alpha1*N1 + D1)*(1-eta1)
      D2new <- (1-g2)*(alpha2*N2 + D2)*(1-eta2)
      D3new <- (1-g3)*(alpha3*N3 + D3)*(1-eta3)
      D4new <- (1-g4)*(alpha4*N4 + D4)*(1-eta4)
      N1new <- g1*(alpha1*N1 + D1)*(1-eta1)
      N2new <- g2*(alpha2*N2 + D2)*(1-eta2)
      N3new <- g3*(alpha3*N3 + D3)*(1-eta3)
      N4new <- g4*(alpha4*N4 + D4)*(1-eta4)
      Rnew  <- Rvector[t]
      return(c(D1new, D2new, D3new, D4new, N1new, N2new, N3new, N4new, Rnew))
    })
  }
  
  ##  Resource uptake function (Hill function)
  uptake_R <- function(r, R, a, b) {
    return((r*R^a) / (b^a + R^a))
  }
  
  
  
  ####
  #### Simulate model -----------------------------------------------------
  ####
  days <- c(1:days_to_track)
  
  # Get "germination" fractions for each year
  num_spp <- length(parms$r)
  
  
  
  # Loop over seasons
  nmsDNR <- names(DNR)
  dormants <- grep("D", names(DNR))
  NR <- DNR[-dormants] 
  nmsNR <- names(NR)
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
                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, alpha4 = alpha4,
                      eta1 = eta1, eta2 = eta2, eta3 = eta3, eta4 = eta4)
    names(DNR) <- nmsDNR
    NR <- DNR[-dormants]
    names(NR) <- nmsNR
  } # next season
  
  return(saved_outs)
  
} #end simulation function



####
####  SIMULATE THE MODEL --------------------------------------------
####
set.seed(123456789) # Set seed to reproduce random results

## Define vectors of parameters to vary
# Initial conditions

DNR <- rbind(c(D=c(25,0.1,0.1,0.1),N=c(25,0.1,0.1,0.1),R=20)) # open community


sig_e_vec <- c(rep(0,100),seq(0,2.5,length.out = 3000))
timesim <- length(sig_e_vec)
getG <- function(sigE, rho, nTime, num_spp) {
  varcov <- matrix(rep(rho*sigE,num_spp*2), num_spp, num_spp)
  diag(varcov) <- sigE
  if(sigE > 0) { varcov <- Matrix::nearPD(varcov)$mat } # coerce matrix to positive definite
  varcov <- as.matrix(varcov)
  e <- rmvnorm(n = nTime, mean = rep(0,num_spp), sigma = varcov)
  g <- exp(e) / (1+exp(e))
  return(g)
}
rho_sim <- 0
gVec <- matrix(nrow=timesim,ncol=4)
for(i in 1:timesim){
  gVec[i,] <- getG(sigE = sig_e_vec[i], rho = rho_sim, nTime = 1, num_spp = 4)
}



##  Define constant parameters in list
constant_parameters <- list (
  seasons = timesim,               # number of seasons to simulate
  days_to_track = 100,             # number of days to simulate in odSolve
  Rmu = 3,                         # mean resource pulse (on log scale)
  Rsd_annual = 0.0,                # std dev of resource pulses (on log scale)
  # sigE = 0,                      # environmental cue variance
  rho = rho_sim,                   # environmental cue correlation between species
  alpha1 = 0.50,                   # live-to-dormant biomass fraction; spp1
  alpha2 = 0.494,                  # live-to-dormant biomass fraction; spp2
  alpha3 = 0.49,                   # live-to-dormant biomass fraction; spp3
  alpha4 = 0.483,                  # live-to-dormant biomass fraction; spp4
  eta1 = 0.1,                      # dormant mortality; spp1
  eta2 = 0.1,                      # dormant mortality; spp2
  eta3 = 0.1,                      # dormant mortality; spp3
  eta4 = 0.1                       # dormant mortality; spp4
)

grow_parameters <- list (
  r   = c(0.2,0.2,0.2,0.2),   # max growth rate for each species
  a   = c(2,2,2,2),           # rate parameter for Hill function 
  b   = c(2.5,2.5,2.5,2.5),   # shape parameter for Hill function
  eps = c(0.5,0.5,0.5,0.5)    # resource-to-biomass efficiency
)

# Make on long vector of named parameters
constant_param_vec <- c(unlist(constant_parameters), unlist(grow_parameters))
parameters <- c(constant_param_vec, DNR)
dnr_slots <- which(names(parameters) == "")
names(parameters)[dnr_slots] <- colnames(DNR)
outs <- do.call(simulate_model,c(list(gVec=gVec),parameters)) 



####
####  PLOTS FOR INFOGRAPHIC
####
library(ggplot2)
library(reshape2)
library(viridis)
library(ggthemes)
library(zoo)

dosim <- as.data.frame(outs[101:3000,5:8])
colnames(dosim) <- c("Species 1","Species 2","Species 3","Species 4")
dosim$total <- rowSums(dosim)
colnames(dosim)[5] <- "Total Biomass"
dosim$iteration <- 1:2900
dosim_long <- melt(dosim, id.vars = "iteration")

mycv <- function(x) {sd(x) / mean(x)}
rolling_cv <- data.frame(cv = rollapply(dosim$`Total Biomass`, width=100, FUN=mycv, fill=NA),
                         iteration = 1:2900)
ggplot(rolling_cv, aes(x=iteration, y=cv))+
  geom_line()+
  ylab("Rolling CV")+
  xlab("Time")+
  scale_y_continuous(limits=c(0,0.3))+
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        panel.background = element_rect("white"),
        axis.line.x = element_line("grey50"),
        axis.line.y = element_line("grey50"),
        axis.title = element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))
ggsave(paste0("../manuscript/components/open_community_infographic_cv.pdf"), width = 3, height = 3, units = "in")

# mycols <- c("#890584","#b1cb3d","#12205b","#0085dc", rgb(52, 0, 66, maxColorValue = 255))
mycols <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
mycols <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", NA)
# mycols <- c("orange","blue","red", "green", "grey25")
ggplot(dosim_long, aes(x=iteration, y=value, color=variable))+
  geom_line(size=0.35)+
  scale_color_manual(values=mycols)+
  # scale_color_viridis(discrete=TRUE, direction = -1, end = 0.9, name="")+
  # scale_alpha_manual(values = c(0.7,0.7,0.7,0.7,1))+
  xlab("Time")+
  ylab("Biomass")+
  scale_y_continuous(limits=c(0,30))+
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        panel.background = element_rect("white"),
        axis.line.x = element_line("white"),
        axis.line.y = element_line("white"),
        axis.title = element_text(size=16))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  guides(color=FALSE, alpha=FALSE)
ggsave(paste0("../manuscript/components/open_community_infographic.pdf"), width = 8, height = 3, units = "in")




####
####  CLOSED COMMUNITY SIMULATION
####
DNR <- rbind(c(D=c(1,0,0,0),N=c(1,0,0,0),R=20)) # closed community
parameters <- c(constant_param_vec, DNR)
dnr_slots <- which(names(parameters) == "")
names(parameters)[dnr_slots] <- colnames(DNR)
outs_close <- do.call(simulate_model,c(list(gVec=gVec),parameters)) 
dosim_close <- as.data.frame(outs_close[101:3000,5:8])
colnames(dosim_close) <- c("Species 1","Species 2","Species 3","Species 4")
dosim_close$total <- rowSums(dosim_close)
colnames(dosim_close)[5] <- "Total Biomass"
dosim_close$iteration <- 1:2900
dosim_long_close <- melt(dosim_close, id.vars = "iteration")

ggplot(dosim_long_close, aes(x=iteration, y=value, color=variable))+
  geom_line()+
  scale_color_manual(values=mycols)+
  # scale_color_viridis(discrete=TRUE, direction = -1, end = 0.9, name="")+
  # scale_alpha_manual(values = c(0.5,0.5,0.5,0.5,1))+
  xlab("Time")+
  ylab("Biomass")+
  scale_y_continuous(limits=c(0,30))+
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        panel.background = element_rect("white"),
        axis.line.x = element_line("white"),
        axis.line.y = element_line("white"),
        axis.title = element_text(size=16))+
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12))+
  guides(color=FALSE, alpha=FALSE)
ggsave(paste0("../manuscript/components/closed_community_infographic.pdf"), width = 8, height = 3, units = "in")

rolling_cv_close <- data.frame(cv = rollapply(dosim_close$`Total Biomass`, width=100, FUN=mycv, fill=NA),
                               iteration = 1:2900)

ggplot(rolling_cv_close, aes(x=iteration, y=cv))+
  geom_line()+
  ylab("Rolling CV")+
  xlab("Time")+
  scale_y_continuous(limits=c(0,0.3))+
  theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
        panel.background = element_rect("white"),
        axis.line.x = element_line("grey50"),
        axis.line.y = element_line("grey50"),
        axis.title = element_text(size=16))+
  theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))
ggsave(paste0("../manuscript/components/closed_community_infographic_cv.pdf"), width = 3, height = 3, units = "in")


# sig_df <- data.frame(variance = sig_e_vec[101:3000],
#                      years = 101:3000)
# ggplot(sig_df, aes(x=years, y=variance))+
#   geom_point(size=2)+
#   ylab(expression(paste("Environmental variance (", sigma[E]^2, ")")))+
#   xlab("Time")+
#   theme(panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
#         panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
#         panel.background = element_rect("white"),
#         axis.line.x = element_line("grey50"),
#         axis.line.y = element_line("grey50"),
#         axis.title = element_text(size=16))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# ggsave(paste0("../manuscript/components/infographic_variance.png"), width = 8, height = 2, units = "in", dpi = 72)
