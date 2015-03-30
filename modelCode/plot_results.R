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


