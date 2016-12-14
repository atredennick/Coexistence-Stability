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
mycols <- brewer.pal(3, "Set1")
mycols <- viridis(2, begin=0.25, end=0.7)
mycols <- c("#15E7A0", "#13CFE8")
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


####
####  SPECIES RICHNESS - ENVIRONMENTAL VARIABILITY RELATIONSHIP; STORAGE EFFECT
####
### Recreate parameter grid
## Define vectors of parameters to vary
n_sig_e <- 11 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
n_rho <- 11 # Number of seasonal standard deviation levels
rho_vec <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vector

##  Create matrix with all possible combinations of varying parameters
varvars <- expand.grid(sig_e_vec, rho_vec )
names(varvars) <- c("sigE", "rho")

##  Read in simulation results
rho0_storage_effect <- readRDS(paste0(path2results,"storage_effect_4species_rho0.RDS"))
n_sig_e <- 100 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 3, length.out=n_sig_e), n_sig_e) # Make a pretty vector
save_multispp_rho0 <- list() # empty storage list
for(i in 1:length(rho0_storage_effect)){
  tmp <- as.data.frame(rho0_storage_effect[[i]])
  names(tmp) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(rho=0,
                        sigE=sig_e_vec[i],
                        cv=tmp_cv,
                        spprich=tmp_spprich)
  
  save_multispp_rho0 <- rbind(save_multispp_rho0, tmp_out)
}





ggplot(save_multispp_rho0, aes(x=spprich, y=cv))+
  geom_point(shape=21, color="black", fill=mycols[1], size=2)+
  stat_smooth(method="lm", color=mycols[1], se=FALSE, size=0.6)+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  theme_bw()+
  my_theme
fig.width <- 133/2
ggsave(paste0(path2figs,"diversity_stability_relationship_storage_effect.png"), width = fig.width, height = 60, units = "mm", dpi = 600)

ggplot(save_multispp_rho0, aes(x=sigE, y=spprich))+
  geom_point(shape=1, col=mycols[1])+
  ylab("S")+
  xlab(expression(sigma[E]^2))+
  my_theme+
  theme(axis.title.y = element_text(angle=0, face="italic"),
        axis.title.x = element_text(face="italic"),
        panel.border = element_rect(fill=NA),
        plot.background = element_rect(fill = "#EFEFEF"))
ggsave(paste0(path2figs,"diversity_envvar_relationship_storage_effect.png"), width = 60, height = 35, units = "mm", dpi = 600)


####
####  SPECIES RICHNESS - ENVIRONMENTAL VARIABILITY RELATIONSHIP; RELATIVE NONLINEARITY
####
### Recreate parameter grid
## Define vectors of parameters to vary
n_rsd <- 50 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
rsd_vec <- as.data.frame(rsd_vec)
names(rsd_vec) <- "Rsd_annual"

##  Read in simulation results
multispp_rln <- readRDS(paste0(path2results,"relative_nonlinearity_4species.RDS"))
save_multispp <- list() # empty storage list
for(i in 1:length(multispp_rln)){
  tmp <- as.data.frame(multispp_rln[[i]])
  names(tmp) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(rsd=rsd_vec[i,"Rsd_annual"],
                        cv=tmp_cv,
                        spprich=tmp_spprich)
  
  save_multispp <- rbind(save_multispp, tmp_out)
}

ggplot(subset(save_multispp, spprich>0), aes(x=spprich, y=cv))+
  geom_point(shape=21, color="black", fill=mycols[2], size=2)+
  stat_smooth(method="lm", color=mycols[2], se=FALSE, size=0.6)+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  theme_bw()+
  my_theme
fig.width <- 133/2
ggsave(paste0(path2figs,"diversity_stability_relationship_relnonlin.png"), width = fig.width, height = 60, units = "mm", dpi = 600)

ggplot(subset(save_multispp, spprich>0), aes(x=rsd, y=spprich))+
  geom_point(shape=1, col=mycols[2])+
  ylab("S")+
  xlab(expression(sigma[R]^phantom(0)))+
  my_theme+
  theme(axis.title.y = element_text(angle=0, face="italic"),
        axis.title.x = element_text(face="italic"),
        panel.border = element_rect(fill=NA),
        plot.background = element_rect(fill = "#EFEFEF"))
ggsave(paste0(path2figs,"diversity_envvar_relationship_relnonlin.png"), width = 60, height = 35, units = "mm", dpi = 600)


####
####  PRE-DEFINED SPECIES COEXISTENCE, DIVERSITY-STABILITY RELATIONSHIPS
####
##  Storage Effect
##  Read in simulation results
sppco_strg <- readRDS(paste0(path2results,"storageeffect_4species_divstability.RDS"))
save_multispp <- list() # empty storage list
for(i in 1:length(sppco_strg)){
  tmp <- as.data.frame(sppco_strg[[i]])
  names(tmp) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(cv=tmp_cv,
                        spprich=tmp_spprich)
  
  save_multispp <- rbind(save_multispp, tmp_out)
}

ggplot(save_multispp, aes(x=spprich, y=cv))+
  geom_point(shape=21, color="black", fill=mycols[1], size=2)+
  stat_smooth(method="lm", color=mycols[1], se=FALSE, size=0.6)+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  theme_bw()+
  my_theme
ggsave(paste0(path2figs,"diversity_stability_relationship_storage_coexist.png"), width = fig.width, height = 60, units = "mm", dpi = 600)

##  Relative Nonlinearity
##  Read in simulation results
sppco_relnonlin <- readRDS(paste0(path2results,"relnonlin_4species_divstability.RDS"))
save_multispp <- list() # empty storage list
for(i in 1:length(sppco_relnonlin)){
  tmp <- as.data.frame(sppco_relnonlin[[i]])
  names(tmp) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(cv=tmp_cv,
                        spprich=tmp_spprich)
  
  save_multispp <- rbind(save_multispp, tmp_out)
}

save_multispp <- subset(save_multispp, spprich > 0)
ggplot(save_multispp, aes(x=spprich, y=cv))+
  geom_point(shape=21, color="black", fill=mycols[2], size=2)+
  stat_smooth(method="lm", color=mycols[2], se=FALSE, size=0.6)+
  xlab("Number of Species")+
  ylab("Variability of Total\nCommunity Biomass (CV)")+
  theme_bw()+
  my_theme
ggsave(paste0(path2figs,"diversity_stability_relationship_relnonlin_coexist.png"),  width = fig.width, height = 60, units = "mm", dpi = 600)

