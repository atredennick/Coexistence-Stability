rm(list=ls(all.names = TRUE))

source("simulate_model_function_multispecies.R")


rs = c(1,1,1,1)
as = c(2,2,2,2)                    # rate parameter for Hill function 
bs = c(2.5,2.5,2.5,2.5)                  # shape parameter for Hill function
etas = c(0.1, 0.12, 0.14, 0.16)

num_spp <- length(rs)

seasons <- 100
Rmu <- 3                         # mean resource pulse (on log scale)
Rsd_annual <- 0                 # std dev of resource pulses (on log scale)
sigE <- 2                     # environmental cue variance
rho = -1
alphas = rep(0.5, num_spp)
betas = rep(0, num_spp)
thetas = rep(0, num_spp)
epss = rep(0.5, num_spp)
Ds=rep(1,num_spp)
Ns=rep(1, num_spp)
R=20 # initial conditions
nu <- 0

varcov <- matrix(rep(rho*sigE,num_spp*2), num_spp, num_spp)
diag(varcov) <- sigE

if(sigE > 0) { varcov <- Matrix::nearPD(varcov)$mat }


test <- simulate_model(seasons = seasons, days_to_track = 20, Rmu=Rmu, 
                       Rsd_annual = Rsd_annual, sigE = sigE, rho = rho, 
                       alphas = alphas, betas = betas, etas = etas, 
                       thetas = thetas, rs=rs, as = as, bs=bs, 
                       epss=epss, nu=nu, Ds=Ds, Ns=Ns, R=R, varcov)

matplot(test[,(num_spp+1):(num_spp*2)], type="l")
