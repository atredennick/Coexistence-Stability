##  plot_sige_gradient_two_rmus.R: script to collate simulation time series from
##  model runs at two levels of mean resource, 30 levels of environmental 
##  variance (up to sigE=10).
##  This is for the storage effect.
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
results_path       <- "../../simulationResults/SI_results/sige_gradient_two_rmus/"
figures_path       <- "../../manuscript/components/"
seasons_to_exclude <- 500
number_of_files    <- length(grep("*.RDS", list.files(results_path)))



####
####  MY PLOTTING THEME ----
####
my_theme <- theme_few()+
  theme(axis.text          = element_text(size=10, color="grey35"),
        axis.title         = element_text(size=12),
        strip.text         = element_text(size=10, color="grey35"),
        legend.title       = element_text(size=10),
        legend.text        = element_text(size=8, color="grey35"),
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
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20))

n_sig_e              <- 30 # Number of cue variance levels
sig_e_vec            <- pretty(seq(0.1, 10, length.out=n_sig_e), n_sig_e) # Make a pretty vector
rho                  <- as.matrix(c((-1/3),0))
prm                  <- expand.grid(as.matrix(sig_e_vec), rho, c(3,6))
DNR_repped           <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
colnames(prm)        <- c("sigE", "rho", "Rmu")
colnames(DNR_repped) <- colnames(DNR)
prm_full             <- cbind(prm, DNR_repped)

if(dim(prm_full)[1] != number_of_files) { stop("WRONG DIMENSIONS; 
                                               CHECK PARAMETER MATRIX") }



####
####  SUMMARIZE RESULTS ----
####
out_files <- as.data.frame(list.files(results_path)[grep("*.RDS", list.files(results_path))]) 
colnames(out_files) <- "file"
suppressWarnings(
  out_files_sep <- out_files %>%
    separate(file, into = c("sc1", "sc2", "sc3", "sc4", "sc5", "sc6", "id"))
)
out_files <- data.frame(out_files, id = as.numeric(out_files_sep[,"id"])) %>%
  arrange(id)

  
sims_summary <- list() # empty storage
for(i in 1:nrow(out_files)) {
  tmp_sim <- readRDS(paste0(results_path, out_files[i,1]))
  tmp_cv  <- get_cv(tmp_sim, seasons_to_exclude)
  tmp_out <- data.frame(tmp_cv,
                        Rmu  = prm_full[i,"Rmu"],
                        rho  = round(prm_full[i,"rho"],2),
                        sigE = prm_full[i,"sigE"])
  sims_summary <- rbind(sims_summary, tmp_out)
}
if(nrow(sims_summary) != number_of_files) { stop("WRONG DIMENSIONS;
                                                 CHECK OUTPUT") }

i=44
tmp_sim <- readRDS(paste0(results_path, out_files[i,1]))
matplot(tmp_sim[,5:8], type="l")

# Calculate mean CV at different realized richness
cv_means <- sims_summary %>%
  group_by(sigE, rho, num_species, Rmu) %>%
  summarise(avg_cv = mean(cv)) %>%
  filter(num_species == 4)

cv_means$Rmu <- paste0("Rmu = ", cv_means$Rmu)




####
####  PLOT RESULTS ----
####
ggplot()+
  geom_line(data=cv_means, aes(x=sigE, y=avg_cv, color=as.factor(rho)), size=0.7, alpha=0.7)+
  geom_point(data=cv_means, aes(x=sigE, y=avg_cv, color=as.factor(rho)), size=1.5, alpha=1)+
  geom_point(data=cv_means, aes(x=sigE, y=avg_cv, color=as.factor(rho)), size=1.5, alpha=1, shape=1, color="grey40")+
  scale_color_brewer(palette = "Set2", name = expression(rho))+
  # scale_color_viridis(end=0.8, discrete=TRUE, name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of Community Biomass")+
  facet_grid(.~Rmu, scales = "free_y")+
  my_theme+
  theme(legend.position = c(0.93, 0.4))+
  theme(legend.background = element_rect(colour = "lightgrey", size=0.25),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=1)))
ggsave(paste0(figures_path, "SI_storage_effect_two_rmus_fourSpeciesOnly.png"), width = 133, height=90, units = "mm", dpi = 120)



####
####  DIG INTO MEAN AND SDS ----
####
sims_four <- filter(sims_summary, num_species==4)
ggplot(sims_four, aes(x=sigE, y=temporal_mean, color=as.factor(rho)))+
  geom_line()+
  facet_wrap("Rmu", scales = "free_y")+
  my_theme
ggplot(sims_four, aes(x=sigE, y=temporal_stddev, color=as.factor(rho)))+
  geom_line()+
  facet_wrap("Rmu", scales = "free_y")+
  my_theme

library(mvtnorm)
getG <- function(sigE, rho, nTime, num_spp=4) {
  varcov       <- matrix(rep(rho*sigE,num_spp*2), num_spp, num_spp)
  diag(varcov) <- sigE
  if(sigE > 0) { varcov <- Matrix::nearPD(varcov)$mat } # crank through nearPD to fix rounding errors 
  varcov <- as.matrix(varcov)
  e      <- rmvnorm(n = nTime, mean = rep(0,num_spp), sigma = varcov)
  g      <- exp(e) / (1+exp(e))
  return(g)
}

rho <- -(1/3)
nTime <- 5000
sig_e_vec <- seq(0,10,by=0.2)
gs <- matrix(ncol=length(sig_e_vec),nrow=nTime*4)
zeros_ones <- numeric(length(sig_e_vec))
for(i in 1:length(sig_e_vec)){
  gs[,i] <- as.vector(getG(sigE=sig_e_vec[i],rho,nTime))
  zeros_ones[i] <- length(which(gs[,i] < 0.01))
}

boxplot(gs, use.cols = T, xaxt="n")
axis(1, at = c(1:length(sig_e_vec)), labels = sig_e_vec)
plot(sig_e_vec,zeros_ones,ylab="Num. of germination events less than 0.01")

extrema_df <- data.frame(sigE = sig_e_vec, zeros = zeros_ones)
ggplot(extrema_df, aes(x = sigE, y=zeros/5000))+
  geom_col()+
  ylab("Proportion of germination events < 0.01")+
  ggtitle(bquote(rho~~"= -1/3"))+
  my_theme
ggsave("../../manuscript/components/SI_big_sigE_gradient_extreme_events.png", height = 3, width=4, units = "in",dpi = 120)


