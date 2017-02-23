##  plot_richness_variability.R



####
####  LOAD LIBRARIES
####
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plyr)
library(reshape2)
library(synchrony)
library(RColorBrewer)
library(viridis)



####
#### INITIALIZATIONS 
####
# Select path to the results and figures
path2results <- "../simulationResults/"
path2figs <- "../manuscript/components/"

seasons_to_exclude <- 500
mycols <- brewer.pal(3, "Set2")
# mycols <- viridis(2, begin=0.25, end=0.7)
# mycols <- c("#15E7A0", "#13CFE8")
my_theme <- theme(legend.title=element_text(size=8, face="bold"),
                  legend.text=element_text(size=8),
                  legend.background = element_rect(colour = "grey45", size=0.5),
                  legend.key = element_blank(),
                  legend.key.size = unit(0.3, "cm"),
                  panel.grid.major = element_line(colour = "white", linetype = "dotted"),
                  panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
                  strip.background = element_blank())

my_theme <- theme_bw()+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="white"),
        panel.background   = element_rect(fill = "#EFEFEF"),
        axis.text          = element_text(size=10, color="grey35", family = "Arial Narrow"),
        axis.title         = element_text(size=12, family = "Arial Narrow", face = "bold"),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color="black"),
        axis.line.y        = element_line(color="black"),
        strip.background   = element_blank(),
        strip.text         = element_text(size=8, color="grey35", family = "Arial Narrow"),
        legend.title       = element_text(size=8, family = "Arial Narrow"),
        legend.text        = element_text(size=6, color="grey35", family = "Arial Narrow"))

my_theme <- theme_few()+
  theme(axis.text          = element_text(size=6, color="grey35"),
        axis.title         = element_text(size=8),
        strip.text         = element_text(size=8, color="grey35"),
        legend.title       = element_text(size=8),
        legend.text        = element_text(size=6, color="grey35"),
        legend.key.size    = unit(0.3, "cm"))

####
####  SPECIES RICHNESS - ENVIRONMENTAL VARIABILITY RELATIONSHIP; STORAGE EFFECT
####
### Recreate parameter grid
# ## Define vectors of parameters to vary
# n_sig_e <- 11 # Number of cue variance levels
# sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
# n_rho <- 11 # Number of seasonal standard deviation levels
# rho_vec <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vector
# 
# ##  Create matrix with all possible combinations of varying parameters
# varvars <- expand.grid(sig_e_vec, rho_vec )
# names(varvars) <- c("sigE", "rho")

##  Read in simulation results
rho0_storage_effect <- readRDS(paste0(path2results,"storage_effect_8species_regional.RDS"))
n_sig_e <- 100 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 10, length.out=n_sig_e), n_sig_e) # Make a pretty vector
save_multispp_rho0 <- list() # empty storage list
for(i in 1:length(rho0_storage_effect)){
  tmp <- as.data.frame(rho0_storage_effect[[i]])
  names(tmp) <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", 
                  "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8","R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(rho=(-1/8),
                        sigE=sig_e_vec[i],
                        cv=tmp_cv,
                        spprich=tmp_spprich,
                        sdev=sd(tmp_totbiomass),
                        avg=mean(tmp_totbiomass))
  
  save_multispp_rho0 <- rbind(save_multispp_rho0, tmp_out)
}
# write.csv(save_multispp_rho0, "../derivedSimulationStats/storage_effect_8species_regional_div-stab.csv")

avg_cv_per_rich <- ddply(save_multispp_rho0, .(spprich), summarise,
                         avg_cv = mean(cv))

ggplot()+
  geom_jitter(data=save_multispp_rho0, aes(x=spprich, y=cv), shape=21, color="grey40", fill=mycols[1], size=2, width=0.05, alpha=0.5)+
  geom_point(data = avg_cv_per_rich, aes(x=spprich, y=avg_cv), shape=19, color="grey40", size=4)+
  geom_line(data = avg_cv_per_rich, aes(x=spprich, y=avg_cv), color="grey40")+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  theme_bw()+
  my_theme
fig.width <- 133/2
# ggsave(paste0(path2figs,"regional_diversity_stability_storage_effect_8species.png"), width = fig.width, height = 60, units = "mm", dpi = 600)

