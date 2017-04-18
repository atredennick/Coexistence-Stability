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
path2results <- "../../simulationResults/SI_results/"
path2figs <- "../../manuscript/components/"

seasons_to_exclude <- 500
mycols <- brewer.pal(3, "Set2")

my_theme <- theme_few()+
  theme(axis.text          = element_text(size=12, color="grey35"),
        axis.title         = element_text(size=14),
        strip.text         = element_text(size=12, color="grey35"),
        legend.title       = element_text(size=12),
        legend.text        = element_text(size=10, color="grey35"),
        legend.key.size    = unit(0.3, "cm"))

####
####  SPECIES RICHNESS - ENVIRONMENTAL VARIABILITY RELATIONSHIP; STORAGE EFFECT
####
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

regional <- ggplot()+
  geom_jitter(data=save_multispp_rho0, aes(x=spprich, y=cv), shape=21, color="grey40", fill=mycols[1], size=2, width=0.05, alpha=0.5)+
  # geom_point(data = avg_cv_per_rich, aes(x=spprich, y=avg_cv), shape=19, color="grey40", size=2)+
  # geom_line(data = avg_cv_per_rich, aes(x=spprich, y=avg_cv), color="grey40")+
  geom_smooth(data=save_multispp_rho0, aes(x=spprich, y=cv), method="loess", se=FALSE, color=mycols[1])+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  scale_x_continuous(breaks=c(1:8), labels = c(1:8))+
  theme_bw()+
  my_theme
# ggsave(paste0(path2figs,"regional_diversity_stability_storage_effect_8species.png"), width = fig.width, height = 60, units = "mm", dpi = 200)




##############                                      #################
##############  SECOND PART -- FROM HPC SAVED FILES #################
##############                                      #################


sim_files <- list.files(paste0(path2results,"eightspp_local/"))
sim_files <- sim_files[grep("*.RDS", sim_files)]
eight_spp_local <- list()
for(i in 1:length(sim_files)){
  tmp <- as.data.frame(readRDS(paste0(path2results,"eightspp_local/",sim_files[i])))
  names(tmp) <- c("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", 
                  "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8","R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  tmp_out <- data.frame(cv=tmp_cv,
                        spprich=tmp_spprich)
  eight_spp_local <- rbind(eight_spp_local, tmp_out)
}

local <- ggplot(eight_spp_local, aes(x=spprich, y=cv))+
  geom_jitter(shape=21, color="grey40", fill=mycols[1], size=2, width=0.05, alpha=0.5) +
  geom_smooth(method="loess", se=FALSE, color=mycols[1])+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  scale_x_continuous(breaks=c(1:8), labels = c(1:8))+
  theme_bw()+
  my_theme

png(filename = paste0(path2figs,"SI_storage_effect_eightspp_local_regional.png"), width = 8, height=3, units = "in", res=100)
grid.arrange(regional,local,ncol=2)
dev.off()
