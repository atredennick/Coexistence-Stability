##  Script to configure and plot results from simulations.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Date:   3.27.2015


# clear the workspace
rm(list=ls())


####
####  Load libraries -------------------------------------------
####
library('ggplot2'); library('plyr')
require('gtools') 


####
####  Get results from desktop folder and combine --------------
####
directory <- "/Users/atredenn/Desktop/storage_effect_simulations/"
files_to_open <- grep(".rds",list.files(directory))
file_list <- list.files(directory)[files_to_open]
file_list <- mixedsort(file_list) 
num_files <- length(file_list)
long_data <- data.frame(sim=NA,rho=NA,sig_e=NA,sig_r=NA,sd_n=NA,mu_n=NA,cv_n=NA,buffer=NA)
Rsd_vec <- seq(0,3,by=0.25)   
for(file in 1:num_files){
  file_now <- readRDS(paste(directory,file_list[file],sep=""))
  file_now$sig_r <- Rsd_vec[file]
  long_data <- rbind(long_data, file_now)
}
saveRDS(long_data, "storage_effect_simulation_output.RDS")

####
####  Read in simulation results -------------------------------
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
####  Plot results ---------------------------------------------
####
ggplot(sim_summary, aes(x=sig_r, y=cv_n, color=rho, group=as.factor(rho)))+
  geom_line()+
  facet_wrap("sig_e")

ggplot(sim_summary, aes(x=sig_r, y=buffer, color=rho, group=as.factor(rho)))+
  geom_hline(aes(yintercept=1), linetype=2)+
  geom_line()+
  geom_point()+
  facet_wrap("sig_e")+
  ylab(expression(CV[community]/CV[resource]))+
  xlab(expression(paste("Resource Variability (", sigma[r], "; on log scale)")))+
  guides(color = guide_colorbar(title = expression(rho)))+
  theme_bw()

ggplot(sim_summary, aes(x=rho, y=buffer, color=sig_r, group=as.factor(sig_r)))+
  geom_line()+
  facet_wrap("sig_e")

ggplot(sim_summary, aes(x=rho, y=cv_n, color=sig_r, group=as.factor(sig_r)))+
  geom_line()+
  facet_wrap("sig_e")


####
#### RELATIVE NONLINEARITY
####

####
####  Get results from desktop folder and combine --------------
####
directory <- "/Users/atredenn/Desktop/rel_nonlin/"
files_to_open <- grep(".rds",list.files(directory))
file_list <- list.files(directory)[files_to_open]
file_list <- mixedsort(file_list) 
num_files <- length(file_list)
long_data <- data.frame(sim=NA,sig_r=NA,sd_n=NA,mu_n=NA,cv_n=NA,buffer=NA)
Rsd_vec <- seq(0,2,by=0.15)  
for(file in 1:num_files){
  file_now <- readRDS(paste(directory,file_list[file],sep=""))
  file_now$sig_r <- Rsd_vec[file]
  long_data <- rbind(long_data, file_now)
}
saveRDS(long_data, "rel_nonlin_simulation_output.RDS")

####
####  Read in simulation results -------------------------------
####
sim_results <- readRDS("rel_nonlin_simulation_output.RDS")
rows_to_remove <- which(is.na(sim_results$sim)==TRUE)
sim_results <- sim_results[-rows_to_remove,]
sim_summary <- ddply(sim_results, .(sig_r), summarise,
                     sd_n = mean(sd_n),
                     mu_n = mean(mu_n),
                     cv_n = mean(cv_n),
                     buffer = mean(buffer))


####
####  Plot results ---------------------------------------------
####
ggplot(sim_summary, aes(x=sig_r, y=buffer))+
  geom_hline(aes(yintercept=1), linetype=2)+
  geom_line()+
  geom_point()+
  ylab(expression(CV[community]/CV[resource]))+
  xlab(expression(paste("Resource Variability (", sigma[r], "; on log scale)")))+
  guides(color = guide_colorbar(title = expression(rho)))+
  theme_bw()

