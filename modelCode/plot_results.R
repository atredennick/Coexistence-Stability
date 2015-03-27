##  Script to configure and plot results from simulations.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Date:   3.27.2015


# clear the workspace
rm(list=ls())


####
####  Load libraries ------------------------------
####
library('ggplot2'); library('plyr')


####
####  Read in simulation results ------------------
####
sim_results <- readRDS("storage_effect_simulation_output.RDS")
rows_to_remove <- which(is.na(sim_results$sim)==TRUE)
sim_results <- sim_results[-rows_to_remove,]
sim_summary <- ddply(sim_results, .(rho, sig_e, sig_r), summarise,
                     sd_n = mean(sd_n),
                     mu_n = mean(mu_n),
                     cv_n = mean(cv_n),
                     buffer = mean(buffer))

####
####  Plot results --------------------------------
####
ggplot(sim_summary, aes(x=sig_r, y=cv_n, color=rho, group=as.factor(rho)))+
  geom_line()+
  facet_wrap("sig_e")


