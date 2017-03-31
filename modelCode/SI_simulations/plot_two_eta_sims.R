##  plot_two_eta_sims.R: script to collate simulation time series from
##  model runs at two levels of species' correlations, two levels of 
##  dormant mortality rates, and 15 levels of environmental variance.
##  This is for the storage effect.
##
##  Author:       Andrew Tredennick (atredenn@gmail.com)
##  Date created: March 7, 2017

# Clear everything
rm(list=ls(all.names = TRUE))



####
####  PACKAGES, PATH NAMES, ETC. ----
####
library(tidyverse)
library(ggthemes)
library(dplyr)
results_path       <- "../../simulationResults/SI_results/vary_etas_results/"
figures_path       <- "../../manuscript/components/"
seasons_to_exclude <- 500
number_of_files    <- length(grep("*.RDS", list.files(results_path)))



####
####  MY PLOTTING THEME ----
####
my_theme <- theme_few()+
  theme(axis.text          = element_text(size=12, color="grey35"),
        axis.title         = element_text(size=14),
        strip.text         = element_text(size=12, color="grey35"),
        legend.title       = element_text(size=12),
        legend.text        = element_text(size=10, color="grey35"),
        legend.key.size    = unit(0.3, "cm"))



####
####  Function for calculating CV ----
####
get_cv <- function(x, seasons_to_exclude){
  x <- as.data.frame(x)
  names(x)       <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates     <- grep("N", colnames(x))
  tmp_totbiomass <- rowSums(x[seasons_to_exclude:nrow(x),livestates])
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
####  RECREATE PARAMETER GRID ----
####
DNR                  <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
                              c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
                              c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
                              c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))
n_sig_e              <- 15 # Number of cue variance levels, approx.
sig_e_vec            <- pretty(seq(0.1, 2, length.out=n_sig_e), n_sig_e) # Make a pretty vector
rho                  <- as.matrix(c((-1/3),0))
prm                  <- expand.grid(as.matrix(sig_e_vec), rho, 1:4)
DNR_repped           <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
colnames(prm)        <- c("sigE", "rho", "dnr_id")
colnames(DNR_repped) <- colnames(DNR)
prm                  <- cbind(prm, DNR_repped)
prm                  <- subset(prm, select = -c(dnr_id))

loweta                <- 0.2
hieta                 <- 0.8
etas                  <- rbind(c(eta1=loweta, eta2=loweta, eta3=loweta, eta4=loweta),
                               c(eta1=hieta, eta2=hieta, eta3=hieta, eta4=hieta))
prm_repped            <- do.call("rbind", replicate(2, prm,  simplify = FALSE))
etas_repped           <- matrix(rep(etas, each=nrow(prm)), ncol=ncol(etas))
colnames(etas_repped) <- colnames(etas)
prm_full              <- cbind(prm_repped, etas_repped)

if(dim(prm_full)[1] != number_of_files) { stop("WRONG DIMENSIONS; 
                                               CHECK PARAMETER MATRIX") }



####
####  SUMMARIZE RESULTS ----
####
out_files <- as.data.frame(list.files(results_path)[grep("*.RDS", list.files(results_path))]) 
colnames(out_files) <- "file"
suppressWarnings(
  out_files_sep <- out_files %>%
    separate(file, into = c("sc1", "sc2", "sc3", "sc4", "sc5", "id"))
)
out_files <- data.frame(out_files, id = as.numeric(out_files_sep[,"id"])) %>%
  arrange(id)

  
sims_summary <- list() # empty storage
for(i in 1:nrow(out_files)) {
  tmp_sim <- readRDS(paste0(results_path, out_files[i,1]))
  tmp_cv  <- get_cv(tmp_sim, seasons_to_exclude)
  tmp_out <- data.frame(tmp_cv,
                        eta  = prm_full[i,"eta1"],
                        rho  = round(prm_full[i,"rho"],2),
                        sigE = prm_full[i,"sigE"])
  print(tmp_out)
  sims_summary <- rbind(sims_summary, tmp_out)
}
if(nrow(sims_summary) != number_of_files) { stop("WRONG DIMENSIONS;
                                                 CHECK OUTPUT") }

# Calculate mean CV at different realized richness
cv_means <- sims_summary %>%
  group_by(sigE, rho, num_species, eta) %>%
  summarise(avg_cv = mean(cv))

# Find species first occurences for lines
firstones <- cv_means %>%
  group_by(rho, eta, num_species) %>%
  summarise(sigE = min(sigE),
            cv = avg_cv[which.min(sigE)])



####
####  PLOT RESULTS ----
####
ggplot()+
  geom_segment(data=firstones, aes(x=sigE, y=cv, xend=sigE, yend=0, color=as.factor(num_species)), size=0.5, alpha=0.7)+
  geom_line(data=cv_means, aes(x=sigE, y=avg_cv, color=as.factor(num_species)), size=0.7, alpha=0.7)+
  geom_point(data=cv_means, aes(x=sigE, y=avg_cv, color=as.factor(num_species)), size=1.5, alpha=1)+
  geom_point(data=cv_means, aes(x=sigE, y=avg_cv, color=as.factor(num_species)), size=1.5, alpha=1, shape=1, color="grey40")+
  scale_color_brewer(palette = "Set2", name = "Species\nRichness")+
  # scale_color_viridis(end=0.8, discrete=TRUE, name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of Community Biomass")+
  facet_grid(eta~rho, scales = "free_y")+
  my_theme+
  theme(legend.position = c(0.93, 0.155))+
  theme(legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=6),
        legend.background = element_rect(colour = "lightgrey", size=0.25),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=1)))
ggsave(paste0(figures_path, "SI_storage_effect_two_etas.png"), width = 6, height=4, units = "in", dpi =120)


