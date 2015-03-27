##  Script to configure and plot results from simulations.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Date:   3.27.2015


# clear the workspace
rm(list=ls())


####
####  Load libraries ------------------------------
####
library('ggplot2'); library('dplyr')


####
####  Read in simulation results ------------------
####
sim_results <- readRDS("storage_effect_simulation_output.RDS")
