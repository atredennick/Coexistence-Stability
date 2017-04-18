

# Clear the memory
rm(list=ls())


####
####  LOAD LIBRARIES
####
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(grid)
library(tidyverse)
library(dplyr)

####
####  PATHS, ETC. ----
####
# Select path to the results and figures
results_path <- "../../simulationResults/"
figure_path  <- "../../manuscript/components/"

seasons_to_exclude <- 500



####
####  HELPER FUNCTIONS FOR PLOTTING ----
####
##  Function for multi-plots with single legend
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="top"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

##  Set color scheme function
get_colors <- function(n_cols){
  return(scales::seq_gradient_pal("red", "blue", "Lab")(seq(0,1,length.out=n_cols)))
}



####
####  PLOT STORAGE EFFECT COEXISTENCE STRENGTH ----
####
## Define vectors of parameters to vary
n_sig_e   <- 11 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
n_rho     <- 11 # Number of seasonal standard deviation levels
rho_vec   <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vector

##  Create matrix with all possible combinations of varying parameters
storage_effect_varvars        <- expand.grid(sig_e_vec, rho_vec )
names(storage_effect_varvars) <- c("sigE", "rho")
storage_effect_varvars        <- unique(storage_effect_varvars)

##  Read in equilibrium simulation results and store in data frame
storage_effect_sims <- readRDS(paste0(results_path, "storage_effect_equilibrium_runs.RDS"))
save_seasons        <- list() # empty storage, to promote to dataframe

for(i in 1:length(storage_effect_sims)){
  tmp          <- as.data.frame(storage_effect_sims[[i]])
  names(tmp)   <- c("D1", "D2", "N1", "N2", "R")
  tmp$rho      <- storage_effect_varvars[i,"rho"]
  tmp$sigE     <- storage_effect_varvars[i,"sigE"]
  tmp$simnum   <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out      <- tmp[seasons_to_exclude:nrow(tmp), ]
  save_seasons <- rbind(save_seasons, tmp_out)
}

##  Calculate CV of Total Community Biomass
save_seasons <- save_seasons %>%
  mutate(total_biomass = N1+N2)

strg_eff_comm_cv <- save_seasons %>%
  group_by(rho, sigE, simnum) %>%
  summarise(cv_biomass = sd(total_biomass)/mean(total_biomass))

##  Read in invasion simulation results; calculate invasion growth rates
storage_effect_invasions <- readRDS(paste0(results_path,"storage_effect_invasion_runs.RDS"))
se_invasion_growth_rate  <- list()
rare_abund_inferior <- 10^-10

for(i in 1:length(storage_effect_invasions)){
  tmp      <- storage_effect_invasions[[i]]
  tmpr     <- log(tmp[1:nrow(tmp),2] / rare_abund_inferior)
  meantmpr <- mean(tmpr)
  tmpdf    <- data.frame(rho=storage_effect_varvars[i,"rho"], 
                         sigE=storage_effect_varvars[i,"sigE"], 
                         growth_rate=meantmpr)
  se_invasion_growth_rate <- rbind(se_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
strg_eff_plot_dat <- se_invasion_growth_rate %>%
  left_join(strg_eff_comm_cv) %>%
  filter(sigE > 0) # get rid of single-species communities

##  Make plots
cc_se <- get_colors(length(unique(strg_eff_plot_dat$rho)))
se_cv_plot <- ggplot(strg_eff_plot_dat, aes(x=sigE, y=cv_biomass, color=as.factor(rho)))+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  scale_color_manual(name=expression(paste("Correlation of species' environmental response (",rho,")")), values=cc_se)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

se_gr_plot <- ggplot(strg_eff_plot_dat, aes(x=sigE, y=growth_rate, color=as.factor(rho)))+
  geom_hline(aes(yintercept=0), linetype=2, alpha=0.9)+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("Log invasion growth rate")+
  scale_color_manual(name=expression(paste("Correlation of species' environmental response (",rho,")")), values=cc_se)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

png(paste0(figure_path,"storage_effect_invasion_growth_rate.png"),
    width = 3.5, height=7, units = "in", res=72)
grid_arrange_shared_legend(se_cv_plot, se_gr_plot)
dev.off()



####
####  PLOT RELATIVE NONLINEARITY RESULTS
####
## Define vectors of parameters to vary
n_rsd <- 11 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector

relative_nonlinearity_sims <- readRDS(paste0(results_path,"relative_nonlinearity_equilibrium_runs.RDS"))
save_seasons_rl <- list()
for(i in 1:length(relative_nonlinearity_sims)){
  tmp <- as.data.frame(relative_nonlinearity_sims[[i]])
  names(tmp) <- c("D1", "D2", "N1", "N2", "R")
  tmp$Rsd <- rsd_vec[i]
  tmp$simnum <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out <- tmp[seasons_to_exclude:nrow(tmp), ]
  save_seasons_rl <- rbind(save_seasons_rl, tmp_out)
}

##  Calculate CV of Total Community Biomass
rl_community_biomass_cv <- save_seasons_rl %>%
  mutate(total_biomass = N1 + N2) %>%
  group_by(Rsd, simnum) %>%
  summarise(cv_biomass = sd(total_biomass)/mean(total_biomass))

##  Calculate Invasion Growth Rate
relative_nonlinearity_invasions <- readRDS(paste0(results_path,"relative_nonlinearity_invasion_runs.RDS"))
rl_invasion_growth_rate <- list()
for(i in 1:length(relative_nonlinearity_invasions)){
  tmp <- relative_nonlinearity_invasions[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),2] / rare_abund_inferior)
  meantmpr <- mean(tmpr)
  tmpdf <- data.frame(Rsd=rsd_vec[i],
                      growth_rate=meantmpr)
  rl_invasion_growth_rate <- rbind(rl_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
rl_plot_dat <- merge(rl_invasion_growth_rate, rl_community_biomass_cv)
# rl_plot_dat <- subset(rl_plot_dat, sigE>0 & Rsd!=1.5)

##  Make the plots
cc_rl <- get_colors(length(unique(rl_plot_dat$sigE)))
rl_cv_plot <- ggplot(rl_plot_dat, aes(x=Rsd, y=cv_biomass))+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3, color="grey35") +
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  theme_few()

rl_gr_plot <- ggplot(rl_plot_dat, aes(x=Rsd, y=growth_rate))+
  geom_hline(aes(yintercept=0), linetype=2, alpha=0.9)+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3, color="grey35")+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("Log invasion growth rate")+
  theme_few()

png(paste0(figure_path,"relative_nonlinearity_invasion_growth_rate.png"),
    width = 3.5, height=7, units = "in", res=72)
grid.arrange(rl_cv_plot, rl_gr_plot)
dev.off()


