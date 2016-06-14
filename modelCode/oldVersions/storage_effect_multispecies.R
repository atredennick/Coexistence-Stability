rm(list=ls(all.names = TRUE))

source("simulate_model_function_multispecies.R")


rs = c(1,1,1,1)
as = c(2,2,2,2)                    # rate parameter for Hill function 
bs = c(2.5,2.5,2.5,2.5)                  # shape parameter for Hill function
etas = c(0.1, 0.1, 0.1, 0.1)

num_spp <- length(rs)

seasons <- 500
Rmu <- 3                         # mean resource pulse (on log scale)
Rsd_annual <- 0                 # std dev of resource pulses (on log scale)
sigE <- 0                   # environmental cue variance
rho = 0
alphas = c(0.5, 0.49, 0.48, 0.47)
betas = rep(0, num_spp)
thetas = rep(0, num_spp)
epss = rep(0.5, num_spp)
Ds=rep(1,num_spp)
Ns=rep(1, num_spp)
R=20 # initial conditions
nu <- 0

varcov <- matrix(rep(rho*sigE,num_spp*2), num_spp, num_spp)
diag(varcov) <- sigE

if(sigE > 0) { varcov <- Matrix::nearPD(varcov)$mat } # coerce matrix to positive definite


test <- simulate_model(seasons = seasons, days_to_track = 20, Rmu=Rmu, 
                       Rsd_annual = Rsd_annual, sigE = sigE, rho = rho, 
                       alphas = alphas, betas = betas, etas = etas, 
                       thetas = thetas, rs=rs, as = as, bs=bs, 
                       epss=epss, nu=nu, Ds=Ds, Ns=Ns, R=R, varcov)

comm_abund <- rowSums(test[,(num_spp+1):(num_spp*2)])

par(mfrow=c(1,2), las=1)
matplot(test[,(num_spp+1):(num_spp*2)], type="l", main=paste("Stability = ", round(sd(comm_abund)/mean(comm_abund),2)))
# lines(comm_abund, lwd=2, col="skyblue")
boxplot(test[,(num_spp+1):(num_spp*2)], col=c(1:num_spp))


