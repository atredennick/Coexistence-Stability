##  plot_means_sdevs_div_var.R: plots community level temporal means and
##  standard deviations at each level of richness at a given level of 
##  environmental variability. This is for storage effect AND relative
##  nonlinearity.

##  Clear the workspace,
rm(list=ls(all.names = TRUE))



####
####  LIBRARIES, PATHS, ETC. ----
####
library(tidyverse)
library(dplyr)
library(ggthemes)
library(gridExtra)

results_path <- "../derivedSimulationStats/"
figures_path <- "../manuscript/components/"



####
####  MY PLOTTING THEME ----
####
my_theme <- theme_few()+
  theme(axis.text          = element_text(size=6, color="grey35"),
        axis.title         = element_text(size=8),
        strip.text         = element_text(size=8, color="grey35"),
        legend.title       = element_text(size=8),
        legend.text        = element_text(size=6, color="grey35"),
        legend.key.size    = unit(0.3, "cm"),
        plot.title         = element_text(size = 8))



####
####  STORAGE EFFECT ----
####
storage_effect_metrics <- read_csv(paste0(results_path, "storage_effect_local_div-stab.csv"))

g1 <- ggplot(storage_effect_metrics, aes(x=as.factor(spprich), y=avg))+
  geom_boxplot()+
  ylab("Temporal Mean of Total Biomass")+
  xlab("Species Richness") +
  ggtitle("Storage Effect:\nTemporal Mean Biomass")+
  my_theme
g2 <- ggplot(storage_effect_metrics, aes(x=as.factor(spprich), y=sdev))+
  geom_boxplot()+
  ylab("Temporal Std. Dev. of Total Biomass")+
  xlab("Species Richness") +
  ggtitle("Storage Effect:\nTemporal Std. Dev. Biomass")+
  my_theme



####
####  RELATIVE NONLINEARITY ----
####
relative_nonlinearity_metrics <- read_csv(paste0(results_path, "relative_nonlinearity_local_div-stab.csv"))

g3 <- ggplot(relative_nonlinearity_metrics, aes(x=as.factor(spprich), y=avg))+
  geom_boxplot()+
  ylab("Temporal Mean of Total Biomass")+
  xlab("Species Richness") +
  ggtitle("Relative Nonlinearity:\nTemporal Mean Biomass")+
  my_theme
g4 <- ggplot(relative_nonlinearity_metrics, aes(x=as.factor(spprich), y=sdev))+
  geom_boxplot()+
  ylab("Temporal Std. Dev. of Total Biomass")+
  xlab("Species Richness") +
  ggtitle("Relative Nonlinearity:\nTemporal Std. Dev. Biomass")+
  my_theme



####
####  COMBINE PLOTS AND SAVE ----
####
png(paste0(figures_path, "SI_temp_means_sds_local.png"), width = 5, height = 5, units = "in", res = 120)
grid.arrange(g1,g2,g3,g4, ncol=2, nrow=2)
dev.off()

