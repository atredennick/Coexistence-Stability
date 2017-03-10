##  plot_relnonlin_three_rmus.R: script to collate simulation time series from
##  model runs at three levels of mean resource, 11 levels of environmental 
##  variance. This is stable to unstable species addition.
##  This is for the relative nonlinearity.
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: March 8, 2017

# Clear everything
rm(list=ls(all.names = TRUE))



####
####  PACKAGES, PATH NAMES, ETC. ----
####
library(tidyverse)
library(ggthemes)
library(dplyr)
results_path1       <- "../../simulationResults/SI_results/relnonlin_four_rmus/"
results_path2       <- "../../simulationResults/SI_results/relnonlin_four_rmus_unstable2stable/"
figures_path        <- "../../manuscript/components/"
seasons_to_exclude  <- 500
number_of_files1    <- length(grep("*.RDS", list.files(results_path1)))
number_of_files2    <- length(grep("*.RDS", list.files(results_path2)))



####
####  MY PLOTTING THEME ----
####
my_theme <- theme_few()+
  theme(axis.text          = element_text(size=6, color="grey35"),
        axis.title         = element_text(size=8),
        strip.text         = element_text(size=8, color="grey35"),
        legend.title       = element_text(size=8),
        legend.text        = element_text(size=6, color="grey35"),
        legend.key.size    = unit(0.3, "cm"))



####
####  Function for calculating CV ----
####
get_cv <- function(x, seasons_to_exclude){
  x <- as.data.frame(x)
  names(x)       <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates     <- grep("N", colnames(x))
  tmp_totbiomass <- rowSums(x[seasons_to_exclude:nrow(x),livestates], na.rm = TRUE)
  tmp_cv         <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_mean       <- mean(tmp_totbiomass)
  tmp_sd         <- sd(tmp_totbiomass)
  tmp_sppavg     <- colMeans(x[seasons_to_exclude:nrow(x),livestates])
  tmp_spprich    <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(cv              = tmp_cv, 
                        num_species     = tmp_spprich, 
                        temporal_mean   = tmp_mean, 
                        temporal_stddev = tmp_sd)
}



####
####  RECREATE PARAMETER GRID (STABLE TO USTABLE) ----
####
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))

## Define vectors of parameters to vary
n_rsd <- 25 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0.1, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
prm <- expand.grid(as.matrix(rsd_vec), 1:4, 1:4)
DNR_repped <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
colnames(prm) <- c("Rsd_annual", "Rmu", "dnr_id")
colnames(DNR_repped) <- colnames(DNR)
prm <- cbind(prm, DNR_repped)
prm_full <- subset(prm, select = -c(dnr_id))


if(dim(prm_full)[1] != number_of_files1) { stop("WRONG DIMENSIONS; 
                                               CHECK PARAMETER MATRIX") }



####
####  SUMMARIZE RESULTS ----
####
out_files <- as.data.frame(list.files(results_path1)[grep("*.RDS", list.files(results_path1))]) 
colnames(out_files) <- "file"
suppressWarnings(
  out_files_sep <- out_files %>%
    separate(file, into = c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", "sc7", "id"))
)
out_files <- data.frame(out_files, id = as.numeric(out_files_sep[,"id"])) %>%
  arrange(id)

  
sims_summary <- list() # empty storage
for(i in 1:nrow(out_files)) {
  tmp_sim <- readRDS(paste0(results_path1, out_files[i,1]))
  tmp_cv  <- get_cv(tmp_sim, seasons_to_exclude)
  tmp_out <- data.frame(tmp_cv,
                        Rmu  = prm_full[i,"Rmu"],
                        Rsd  = prm_full[i,"Rsd_annual"])
  sims_summary <- rbind(sims_summary, tmp_out)
}
if(nrow(sims_summary) != number_of_files1) { stop("WRONG DIMENSIONS;
                                                 CHECK OUTPUT") }

# Calculate mean CV at different realized richness
torms <- which(is.nan(sims_summary$cv))
sims_summary <- sims_summary[-torms,]
cv_means_su <- sims_summary %>%
  group_by(Rmu, Rsd, num_species) %>%
  summarise(avg_cv = mean(cv))  %>%
  filter(num_species>0 & avg_cv<5 & avg_cv>0)

cv_means_su$Rmu <- paste0("Rmu = ", cv_means_su$Rmu)

# Find species first occurences for lines
firstones_su <- cv_means_su %>%
  group_by(Rmu, num_species) %>%
  summarise(Rsd = min(Rsd),
            cv = avg_cv[which.min(Rsd)])






####
####  RECREATE PARAMETER GRID (USTABLE TO STABLE) ----
####
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))

## Define vectors of parameters to vary
n_rsd <- 25 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0.1, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
prm <- expand.grid(as.matrix(rsd_vec), 1:4, 1:4)
DNR_repped <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
colnames(prm) <- c("Rsd_annual", "Rmu", "dnr_id")
colnames(DNR_repped) <- colnames(DNR)
prm <- cbind(prm, DNR_repped)
prm_full <- subset(prm, select = -c(dnr_id))


if(dim(prm_full)[1] != number_of_files2) { stop("WRONG DIMENSIONS; 
                                               CHECK PARAMETER MATRIX") }



####
####  SUMMARIZE RESULTS ----
####
out_files <- as.data.frame(list.files(results_path2)[grep("*.RDS", list.files(results_path2))]) 
colnames(out_files) <- "file"
suppressWarnings(
  out_files_sep <- out_files %>%
    separate(file, into = c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", "sc7", "id"))
)
out_files <- data.frame(out_files, id = as.numeric(out_files_sep[,"id"])) %>%
  arrange(id)

# which(!(1:432 %in% out_files$id))

sims_summary <- list() # empty storage
for(i in 1:nrow(out_files)) {
  tmp_sim <- readRDS(paste0(results_path2, out_files[i,1]))
  tmp_cv  <- get_cv(tmp_sim, seasons_to_exclude)
  tmp_out <- data.frame(tmp_cv,
                        Rmu  = prm_full[out_files[i,2],"Rmu"],
                        Rsd  = prm_full[out_files[i,2],"Rsd_annual"])
  sims_summary <- rbind(sims_summary, tmp_out)
}
if(nrow(sims_summary) != number_of_files2) { stop("WRONG DIMENSIONS;
                                                 CHECK OUTPUT") }

# Calculate mean CV at different realized richness
cv_means <- sims_summary %>%
  group_by(Rmu, Rsd, num_species) %>%
  summarise(avg_cv = mean(cv)) %>%
  filter(num_species>0 & avg_cv<5 & avg_cv>0)

cv_means$Rmu <- paste0("Rmu = ", cv_means$Rmu)

# Find species first occurences for lines
firstones <- cv_means %>%
  group_by(Rmu, num_species) %>%
  summarise(Rsd = min(Rsd),
            cv = avg_cv[which.min(Rsd)])



####
####  COMBINE RESULTS ----
####
cv_means$spporder <- "unstable2stable"
cv_means_su$spporder <- "stable2unsable"
firstones$spporder <- "unstable2stable"
firstones_su$spporder <- "stable2unsable"

cv_means <- rbind(cv_means, cv_means_su)
firstones <- rbind(firstones, firstones_su)


####
####  PLOT RESULTS ----
####
ggplot()+
  geom_segment(data=firstones, aes(x=Rsd, y=cv, xend=Rsd, yend=0, color=as.factor(num_species)), size=0.5, alpha=0.7)+
  geom_line(data=cv_means, aes(x=Rsd, y=avg_cv, color=as.factor(num_species)), size=0.7, alpha=0.7)+
  geom_point(data=cv_means, aes(x=Rsd, y=avg_cv, color=as.factor(num_species)), size=1.5, alpha=1)+
  geom_point(data=cv_means, aes(x=Rsd, y=avg_cv, color=as.factor(num_species)), size=1.5, alpha=1, shape=1, color="grey40")+
  scale_color_brewer(palette = "Set2", name = "Species\nRichness")+
  xlab(expression(paste("Standard deviation of resource (",sigma[R], ")")))+
  ylab("CV of Community Biomass")+
  facet_grid(spporder~Rmu, scales = "free_y")+
  my_theme+
  # theme(legend.position = c(0.93, 0.155))+
  theme(legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=6),
        legend.background = element_rect(colour = "lightgrey", size=0.25),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=1)))
# ggsave(paste0(figures_path, "SI_relative_nonlinearity_four_rmus.png"), width = 160, height=60, units = "mm", dpi = 600)
