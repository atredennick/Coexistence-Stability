

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

# png(paste0(path2figs,"storage_effect_factorial_plots.png"),
    # width = 3.5, height=7, units = "in", res=72)
grid_arrange_shared_legend(se_cv_plot, se_gr_plot)
# dev.off()


####
####  PLOT STORAGE EFFECT RESULTS -- VARIABLE RESOURCE
####
## Define vectors of parameters to vary
n_rsd <- 11 # Number of cue variance levels
rsd_vec <- pretty(seq(0, 1.4, length.out=n_sig_e), n_rsd) # Make a pretty vector
n_rho <- 11 # Number of seasonal standard deviation levels
rho_vec <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vector
##  Create matrix with all possible combinations of varying parameters
storage_effect_varvars2 <- expand.grid(rsd_vec, rho_vec )
names(storage_effect_varvars2) <- c("rsd", "rho")
storage_effect_varvars2 <- unique(storage_effect_varvars2)

storage_effect_sims2 <- readRDS(paste0(path2results,"storage_effect_equilibrium_varyresource_runs.RDS"))
save_seasons_se2 <- list()
for(i in 1:length(storage_effect_sims2)){
  tmp <- as.data.frame(storage_effect_sims2[[i]])
  names(tmp) <- c("D1", "D2", "N1", "N2", "R")
  tmp$rho <- storage_effect_varvars2[i,"rho"]
  tmp$rsd <- storage_effect_varvars2[i,"rsd"]
  tmp$simnum <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out <- tmp[seasons_to_exclude:nrow(tmp), ]
  save_seasons_se2 <- rbind(save_seasons_se2, tmp_out)
}

##  Calculate CV of Total Community Biomass
save_seasons_se2$total_biomass <- with(save_seasons_se2, N1+N2)
se2_community_biomass_cv <- ddply(save_seasons_se2, .(rho, rsd, simnum), summarise,
                                 cv_biomass = sd(total_biomass)/mean(total_biomass))


##  Calculate Invasion Growth Rate
storage_effect_invasions2 <- readRDS(paste0(path2results,"storage_effect_invasion_varyresource_runs.RDS"))
se2_invasion_growth_rate <- list()
for(i in 1:length(storage_effect_invasions2)){
  tmp <- storage_effect_invasions2[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),2] / 1)
  meantmpr <- mean(tmpr)
  tmpdf <- data.frame(rho=storage_effect_varvars2[i,"rho"], rsd=storage_effect_varvars2[i,"rsd"], growth_rate=meantmpr)
  se2_invasion_growth_rate <- rbind(se2_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
se2_plot_dat <- merge(se2_invasion_growth_rate, se2_community_biomass_cv)

##  Make plots
cc_se2 <- get_colors(length(unique(se2_plot_dat$rho)))
se2_cv_plot <- ggplot(se2_plot_dat, aes(x=rsd, y=cv_biomass, color=as.factor(rho)))+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  scale_color_manual(name=expression(paste("Correlation of species' environmental response (",rho,")")), values=cc_se)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

se2_gr_plot <- ggplot(se2_plot_dat, aes(x=rsd, y=growth_rate, color=as.factor(rho)))+
  geom_hline(aes(yintercept=0), linetype=2, alpha=0.9)+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("Log invasion growth rate")+
  scale_y_continuous(breaks=c(-0.1,0,0.1), limits=c(-0.1,0.13))+
  scale_color_manual(name=expression(paste("Correlation of species' environmental response (",rho,")")), values=cc_se)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

png(paste0(path2figs,"storage_effect_factorial_varyresource_plots.png"),
    width = 3.5, height=7, units = "in", res=72)
grid_arrange_shared_legend(se2_cv_plot, se2_gr_plot)
dev.off()



####
####  PLOT STORAGE EFFECT RESULTS -- VARIABLE RESOURCE, VARIABLE CUE
####
## Define vectors of parameters to vary
n_rsd <- 11 # Number of cue variance levels
rsd_vec <- pretty(seq(0, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
n_sig_e <- 11 # Number of seasonal standard deviation levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
##  Create matrix with all possible combinations of varying parameters
storage_effect_varvars3 <- expand.grid(rsd_vec, sig_e_vec)
names(storage_effect_varvars3) <- c("rsd", "sigE")
storage_effect_varvars3 <- unique(storage_effect_varvars3)

storage_effect_sims3 <- readRDS(paste0(path2results,"storage_effect_equilibrium_varyresource_varycue_runs.RDS"))
save_seasons_se3 <- list()
for(i in 1:length(storage_effect_sims3)){
  tmp <- as.data.frame(storage_effect_sims3[[i]])
  names(tmp) <- c("D1", "D2", "N1", "N2", "R")
  tmp$sigE <- storage_effect_varvars3[i,"sigE"]
  tmp$rsd <- storage_effect_varvars3[i,"rsd"]
  tmp$simnum <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out <- tmp[seasons_to_exclude:nrow(tmp), ]
  save_seasons_se3 <- rbind(save_seasons_se3, tmp_out)
}

##  Calculate CV of Total Community Biomass
save_seasons_se3$total_biomass <- with(save_seasons_se3, N1+N2)
se3_community_biomass_cv <- ddply(save_seasons_se3, .(sigE, rsd, simnum), summarise,
                                  cv_biomass = sd(total_biomass)/mean(total_biomass))


##  Calculate Invasion Growth Rate
storage_effect_invasions3 <- readRDS(paste0(path2results,"storage_effect_invasion_varyresource_varycue_runs.RDS"))
se3_invasion_growth_rate <- list()
for(i in 1:length(storage_effect_invasions3)){
  tmp <- storage_effect_invasions3[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),2] / 1)
  meantmpr <- mean(tmpr)
  tmpdf <- data.frame(sigE=storage_effect_varvars3[i,"sigE"], rsd=storage_effect_varvars3[i,"rsd"], growth_rate=meantmpr)
  se3_invasion_growth_rate <- rbind(se3_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
se3_plot_dat <- merge(se3_invasion_growth_rate, se3_community_biomass_cv)
se3_plot_dat <- subset(se3_plot_dat, cv_biomass<10)

##  Make plots
cc_se3 <- get_colors(length(unique(se3_plot_dat$sigE)))
se3_cv_plot <- ggplot(se3_plot_dat, aes(x=rsd, y=cv_biomass, color=as.factor(sigE)))+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  scale_color_manual(name=expression(paste("Variance of environmental cue (",sigma[E]^2, ")")), values=cc_se3)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

se3_gr_plot <- ggplot(se3_plot_dat, aes(x=rsd, y=growth_rate, color=as.factor(sigE)))+
  geom_hline(aes(yintercept=0), linetype=2, alpha=0.9)+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("Log invasion growth rate")+
  # scale_y_continuous(breaks=c(-0.1,0,0.1), limits=c(-0.1,0.13))+
  scale_color_manual(name=expression(paste("Variance of environmental cue (",sigma[E]^2, ")")), values=cc_se3)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

png(paste0(path2figs,"storage_effect_factorial_varyresource_varycue_plots.png"),
    width = 3.5, height=7, units = "in", res=72)
grid_arrange_shared_legend(se3_cv_plot, se3_gr_plot)
dev.off()




####
####  PLOT RELATIVE NONLINEARITY RESULTS
####
## Define vectors of parameters to vary
n_sig_e <- 11 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
n_rsd <- 11 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0, 1.5, length.out=n_rsd), n_rsd) # Make a pretty vector

##  Create matrix with all possible combinations of varying parameters
relative_nonlinearity_varvars <- expand.grid(sig_e_vec, rsd_vec )
names(relative_nonlinearity_varvars) <- c("sigE", "Rsd")
relative_nonlinearity_varvars <- unique(relative_nonlinearity_varvars)

relative_nonlinearity_sims <- readRDS(paste0(path2results,"relative_nonlinearity_equilibrium_runs.RDS"))
save_seasons_rl <- list()
for(i in 1:length(relative_nonlinearity_sims)){
  tmp <- as.data.frame(relative_nonlinearity_sims[[i]])
  names(tmp) <- c("D1", "D2", "N1", "N2", "R")
  tmp$Rsd <- relative_nonlinearity_varvars[i,"Rsd"]
  tmp$sigE <- relative_nonlinearity_varvars[i,"sigE"]
  tmp$simnum <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out <- tmp[seasons_to_exclude:nrow(tmp), ]
  save_seasons_rl <- rbind(save_seasons_rl, tmp_out)
}

##  Calculate CV of Total Community Biomass
save_seasons_rl$total_biomass <- with(save_seasons_rl, N1+N2)
rl_community_biomass_cv <- ddply(save_seasons_rl, .(Rsd, sigE, simnum), summarise,
                                 cv_biomass = sd(total_biomass)/mean(total_biomass))

##  Calculate Invasion Growth Rate
relative_nonlinearity_invasions <- readRDS(paste0(path2results,"relative_nonlinearity_invasion_runs.RDS"))
rl_invasion_growth_rate <- list()
for(i in 1:length(relative_nonlinearity_invasions)){
  tmp <- relative_nonlinearity_invasions[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),2] / 1)
  meantmpr <- mean(tmpr)
  tmpdf <- data.frame(Rsd=relative_nonlinearity_varvars[i,"Rsd"], sigE=relative_nonlinearity_varvars[i,"sigE"], growth_rate=meantmpr)
  rl_invasion_growth_rate <- rbind(rl_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
rl_plot_dat <- merge(rl_invasion_growth_rate, rl_community_biomass_cv)
rl_plot_dat <- subset(rl_plot_dat, sigE>0 & Rsd!=1.5)

##  Make the plots
cc_rl <- get_colors(length(unique(rl_plot_dat$sigE)))
rl_cv_plot <- ggplot(rl_plot_dat, aes(x=Rsd, y=cv_biomass, color=as.factor(sigE)))+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  scale_color_manual(name=expression(paste("Variance of environmental cue (",sigma[E]^2, ")")), values=cc_rl)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

rl_gr_plot <- ggplot(rl_plot_dat, aes(x=Rsd, y=growth_rate, color=as.factor(sigE)))+
  geom_hline(aes(yintercept=0), linetype=2, alpha=0.9)+
  geom_point(size=2, alpha=0.5)+
  geom_point(size=2, alpha=0.8, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=0.3)+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("Log invasion growth rate")+
  scale_color_manual(name=expression(paste("Variance of environmental cue (",sigma[E]^2, ")")), values=cc_rl)+
  theme_few()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))+
  guides(color=guide_legend(nrow=2,byrow=TRUE, title.position = "top"))

png(paste0(path2figs,"relative_nonlinearity_factorial_plots.png"),
    width = 3.5, height=7, units = "in", res=72)
grid_arrange_shared_legend(rl_cv_plot, rl_gr_plot)
dev.off()



####
####  MAKE STORAGE EFFECT AND RELATIVE NONLINEARITY PROFILE PLOTS
####
## Define vectors of parameters to vary
both_varvars <- matrix(c(-0.5,1,2,2.5,5,5,20,   # relative nonlinearity present
                         -0.5,3,5,20,3,5,20,    # rln absent; medium start
                         0.5,1,2,2.5,5,5,20,   # relative nonlinearity present
                         0.5,3,5,20,3,5,20,    # rln absent; medium start
                         -0.5,1,2,2.5,1,2,2.5,  # rln absent; low start
                         0.5,1,2,2.5,1,2,2.5,  # rln absent; low start
                         -0.5,5,5,20,5,5,20,    # rln absent; high start
                         0.5,5,5,20,5,5,20     # rln absent; high start
                        ), ncol=7, byrow = TRUE)
colnames(both_varvars) <- c("rho", "r1", "a1", "b1", "r2", "a2", "b2")
sim_start_names <- c("NOPE", "MED", 
               "NOPE", "MED",
               "LOW", "LOW", 
               "HIGH", "HIGH")

##  Plot the resource uptake curves
mycols <- c("grey", "skyblue")

mytheme <- theme(axis.text=element_text(size=6),
                 axis.title=element_text(size=8),
                 legend.text=element_text(size=6),
                 legend.key.size=unit(0.3, "cm"))

make_Ruptake_plot <- function(df, extratheme){
  require(ggplot2)
  require(ggthemes)
  
  outplot <- ggplot()+
    geom_line(data=df, aes(x=resource, y=uptake, color=as.factor(species), linetype=as.factor(species)))+
    theme_few()+
    scale_color_manual(values=mycols, name="", labels=c("Species A", "Species B"))+
    guides(size=FALSE, linetype=FALSE)+
    scale_y_continuous(limits=c(0,5))+
    xlab("")+
    ylab("")+
    guides(color=FALSE)+
    extratheme
  
  return(outplot)
}

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, alpha, beta){
  return((r*R^alpha) / (beta^alpha + R^alpha))
}

get_Ruptakes <- function(parms){
  R <- seq(0,100,1)
  out_r <- matrix(ncol=2, nrow=length(R))
  alpha <- parms$alpha
  beta <- parms$beta
  for(i in 1:nrow(out_r)){
    out_r[i,1] <- uptake_R(parms$r[1], R[i], alpha[1], beta[1])
    out_r[i,2] <- uptake_R(parms$r[2], R[i], alpha[2], beta[2])
  }
  uptake <- data.frame(species=rep(c(1:2), each=nrow(out_r)),
                       resource=rep(R, times=2),
                       uptake=c(out_r[,1], out_r[,2]))
  return(uptake)
}


##  RELATIVE NONLINEARITY OPERATING
parms <- list(
  r = c(1,5),         # max growth rate for each species
  alpha = c(2,5),     # rate parameter for Hill function 
  beta = c(2.5,20)    # shape parameter for Hill function
)

relnonlin_uptake <- get_Ruptakes(parms)
rel_on_plot <- make_Ruptake_plot(relnonlin_uptake, mytheme)
ggsave("../manuscript/components/Ruptake_relnonlin.png", plot = rel_on_plot, 
       width = 2, height = 1, units = "in", dpi = 72)

# Start low
parms <- list(
  r = c(1,1),         # max growth rate for each species
  alpha = c(2,2),     # rate parameter for Hill function 
  beta = c(2.5,2.5)    # shape parameter for Hill function
)
low_uptake <- get_Ruptakes(parms)
start_lo_plot <- make_Ruptake_plot(low_uptake, mytheme)
ggsave("../manuscript/components/start_lo_plot.png", plot = start_lo_plot, 
       width = 2, height = 1, units = "in", dpi = 72)

# Start med
parms <- list(
  r = c(3,3),         # max growth rate for each species
  alpha = c(5,5),     # rate parameter for Hill function 
  beta = c(20,20)    # shape parameter for Hill function
)
med_uptake <- get_Ruptakes(parms)
start_med_plot <- make_Ruptake_plot(med_uptake, mytheme)
ggsave("../manuscript/components/start_med_plot.png", plot = start_med_plot, 
       width = 2, height = 1, units = "in", dpi = 72)

# Start high
parms <- list(
  r = c(5,5),         # max growth rate for each species
  alpha = c(5,5),     # rate parameter for Hill function 
  beta = c(20,20)    # shape parameter for Hill function
)
high_uptake <- get_Ruptakes(parms)
start_hi_plot <- make_Ruptake_plot(high_uptake, mytheme)
ggsave("../manuscript/components/start_hi_plot.png", plot = start_hi_plot, 
       width = 2, height = 1, units = "in", dpi = 72)



##  Get simulation results and collate
both_sims <- readRDS(paste0(path2results,"storageeffect_relativenonlinearity_equilibrium_runs.RDS"))
save_seasons_both <- list()
for(i in 1:length(both_sims)){
  tmp <- as.data.frame(both_sims[[i]])
  names(tmp) <- c("D1", "D2", "N1", "N2", "R")
  tmp$rho <- both_varvars[i,"rho"]
  if(both_varvars[i,"a1"] == both_varvars[i,"a2"]) rel_nonlin <- "no"
  if(both_varvars[i,"a1"] != both_varvars[i,"a2"]) rel_nonlin <- "yes"
  tmp$rel_nonlin <- rel_nonlin
  tmp$simnum <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out <- tmp[seasons_to_exclude:nrow(tmp), ]
  tmp_out$sim_code <- sim_start_names[i]
  save_seasons_both <- rbind(save_seasons_both, tmp_out)
}

##  Calculate CV of Total Community Biomass
save_seasons_both$total_biomass <- with(save_seasons_both, N1+N2)
both_community_biomass_cv <- ddply(save_seasons_both, .(rho, rel_nonlin, sim_code), summarise,
                                 cv_biomass = sd(total_biomass)/mean(total_biomass))
tmpdf <- rbind(both_community_biomass_cv[which(both_community_biomass_cv$sim_code=="NOPE"),],
               both_community_biomass_cv[which(both_community_biomass_cv$sim_code=="NOPE"),],
               both_community_biomass_cv[which(both_community_biomass_cv$sim_code=="NOPE"),])
tmpdf$sim_code <- c(rep("LOW",2), rep("MED",2), rep("HIGH",2))
both_community_biomass_cv <- rbind(tmpdf, both_community_biomass_cv[which(both_community_biomass_cv$sim_code!="NOPE"),])

make_cv_plot <- function(df, ylims){
  require(ggplot2)
  require(ggthemes)
  
  ggplot(df, aes(x=as.factor(rho),y=cv_biomass, group=rel_nonlin))+
    geom_line()+
    geom_point(size=4,color="white")+
    geom_point(size=4, aes(shape=rel_nonlin))+
    ylab("CV of community biomass")+
    xlab(expression(rho))+
    scale_y_continuous(limits=ylims)+
    theme_few()+
    scale_shape_manual(values=c(1,19), name="", labels=c("RNL absent", "RNL present"))+
    theme(legend.position=c(0.75,0.15))+
    theme(legend.text=element_text(size=10),
          legend.key.size=unit(0.5, "cm"))+
    # guides(shape = guide_legend(override.aes = list(size=2)))+
    guides(shape=FALSE)+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=15),
          legend.text=element_text(size=6),
          legend.key.size=unit(0.3, "cm"))
}

ylimits <- c(min(both_community_biomass_cv$cv_biomass), max(both_community_biomass_cv$cv_biomass))
low_gr_rnl_strg_cv <- make_cv_plot(subset(both_community_biomass_cv, sim_code=="LOW"), ylimits)
med_gr_rnl_strg_cv <- make_cv_plot(subset(both_community_biomass_cv, sim_code=="MED"), ylimits)
high_gr_rnl_strg_cv <- make_cv_plot(subset(both_community_biomass_cv, sim_code=="HIGH"), ylimits)


##  Calculate Invasion Growth Rate
both_invasions <- readRDS(paste0(path2results,"storageeffect_relativenonlinearity_invasion_runs.RDS"))
both_invasion_growth_rate <- list()
for(i in 1:length(both_invasions)){
  tmp <- both_invasions[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),2] / 1)
  meantmpr <- mean(tmpr)
  if(both_varvars[i,"a1"] == both_varvars[i,"a2"]) rel_nonlin <- "no"
  if(both_varvars[i,"a1"] != both_varvars[i,"a2"]) rel_nonlin <- "yes"
  tmpdf <- data.frame(rho=both_varvars[i,"rho"], 
                      rel_nonlin=rel_nonlin, 
                      growth_rate=meantmpr)
  tmpdf$sim_code <- sim_start_names[i]
  both_invasion_growth_rate <- rbind(both_invasion_growth_rate, tmpdf)
}
tmpdf <- rbind(both_invasion_growth_rate[which(both_invasion_growth_rate$sim_code=="NOPE"),],
               both_invasion_growth_rate[which(both_invasion_growth_rate$sim_code=="NOPE"),],
               both_invasion_growth_rate[which(both_invasion_growth_rate$sim_code=="NOPE"),])
tmpdf$sim_code <- c(rep("LOW",2), rep("MED",2), rep("HIGH",2))
both_invasion_growth_rate <- rbind(tmpdf, both_invasion_growth_rate[which(both_invasion_growth_rate$sim_code!="NOPE"),])

make_invasion_gr_plot <- function(df, ylims){
  require(ggplot2)
  require(ggthemes)
  
  ggplot(df, aes(x=as.factor(rho),y=growth_rate, group=rel_nonlin))+
    geom_hline(aes(yintercept=0))+
    geom_line()+
    geom_point(size=4,color="white")+
    geom_point(size=4, aes(shape=rel_nonlin))+
    ylab("Log invasion growth rate")+
    xlab(expression(rho))+  
    scale_y_continuous(limits=ylims)+
    theme_few()+
    scale_shape_manual(values=c(19,1))+
    guides(shape=FALSE)+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=15),
          legend.text=element_text(size=6),
          legend.key.size=unit(0.3, "cm"))
}

ylimits <- c(min(both_invasion_growth_rate$growth_rate), max(both_invasion_growth_rate$growth_rate))
low_gr_rnl_strg_gr <- make_invasion_gr_plot(subset(both_invasion_growth_rate, sim_code=="LOW"), ylimits)
med_gr_rnl_strg_gr <- make_invasion_gr_plot(subset(both_invasion_growth_rate, sim_code=="MED"), ylimits)
high_gr_rnl_strg_gr <- make_invasion_gr_plot(subset(both_invasion_growth_rate, sim_code=="HIGH"), ylimits)


png(paste0(path2figs,"storage_effect_relative_nonlinearity_anova_plot_low.png"),
    width = 3.5, height=7, units = "in", res=72)
grid.arrange(low_gr_rnl_strg_cv, low_gr_rnl_strg_gr)
dev.off()

png(paste0(path2figs,"storage_effect_relative_nonlinearity_anova_plot_med.png"),
    width = 3.5, height=7, units = "in", res=72)
grid.arrange(med_gr_rnl_strg_cv, med_gr_rnl_strg_gr)
dev.off()

png(paste0(path2figs,"storage_effect_relative_nonlinearity_anova_plot_high.png"),
    width = 3.5, height=7, units = "in", res=72)
grid.arrange(high_gr_rnl_strg_cv, high_gr_rnl_strg_gr)
dev.off()




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
multispp_strg <- readRDS(paste0(path2results,"storage_effect_4species.RDS"))
save_multispp <- list() # empty storage list
for(i in 1:length(multispp_strg)){
  tmp <- as.data.frame(multispp_strg[[i]])
  names(tmp) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(tmp))
  tmp_totbiomass <- rowSums(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(tmp[seasons_to_exclude:nrow(tmp),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- data.frame(rho=varvars[i,"rho"],
                        sigE=varvars[i,"sigE"],
                        cv=tmp_cv,
                        spprich=tmp_spprich)

  save_multispp <- rbind(save_multispp, tmp_out)
}


library("plot3D")
M <- as.matrix(dcast(save_multispp, rho~sigE, value.var = "cv"))
png(paste0(path2figs,"cv_persp_plot.png"), width = 5, height = 5, units="in", res = 120)
persp3D(x = unique(save_multispp$rho),
       y = unique(save_multispp$sigE),
       z = M[,2:ncol(M)],
       xlab = "", bty = "b2",
       ylab = "", zlab = "", clab = "",
       expand = 0.5, d = 2, phi = 20, theta = 30, resfac = 5,
       contour = list(col = "grey", side = c("zmin", "z")),
       colkey = list(side = 1, length = 0.5),
       ticktype="detailed")
dev.off()


Mspp <- as.matrix(dcast(save_multispp, rho~sigE, value.var = "spprich"))
png(paste0(path2figs,"spprich_persp_plot.png"), width = 5, height = 5, units="in", res = 120)
persp3D(x = unique(save_multispp$rho),
        y = unique(save_multispp$sigE),
        z = Mspp[,2:ncol(M)],
        xlab = "", bty = "b2",
        ylab = "", zlab = "", clab = "",
        expand = 0.5, d = 2, phi = 20, theta = 30, resfac = 5,
        contour = list(col = "grey", side = c("zmin", "z")),
        colkey = list(side = 1, length = 0.5),
        ticktype="detailed")
dev.off()

ggplot(subset(save_multispp, rho!=1), aes(x=spprich, y=cv))+
  geom_point(shape=21, color="white", fill="black", size=3)+
  stat_smooth(method="lm", color="black", fill="grey", se=FALSE)+
  facet_wrap("rho", nrow=2)+
  xlab("Number of Species")+
  ylab("Variability of Total Community Biomass (CV)")+
  theme_few()
# ggsave(paste0(path2figs,"diversity_stability_relationship.png"), width = 8, height = 4, units = "in", dpi = 100)


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
  geom_point(shape=21, color="black", fill="dodgerblue", size=3)+
  stat_smooth(method="lm", color="skyblue", fill="grey", se=FALSE)+
  # stat_smooth(data=subset(save_multispp_rho0, rho==0&spprich>1), aes(x=spprich, y=cv), method="lm", color="black", fill="grey", se=FALSE, linetype=2)+
  xlab("Number of Species")+
  ylab("Variability of Total Community Biomass (CV)")+
  theme_few()
ggsave(paste0(path2figs,"diversity_stability_relationship_storage_effect.png"), width = 4, height = 4, units = "in", dpi = 120)

ggplot(save_multispp_rho0, aes(x=sigE, y=spprich))+
  geom_point()+
  ylab("S")+
  xlab(expression(sigma[E]^2))+
  theme_classic()+
  theme(axis.title.y = element_text(angle=0, face="italic", size=22),
        axis.title.x = element_text(face="italic", size=22))
ggsave(paste0(path2figs,"diversity_envvar_relationship_storage_effect.png"), width = 3, height = 1.5, units = "in", dpi = 72)


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
  geom_point(shape=21, color="black", fill="coral", size=3)+
  stat_smooth(method="lm", color="orange", fill="grey", se=FALSE)+
  xlab("Number of Species")+
  ylab("Variability of Total Community Biomass (CV)")+
  theme_few()
ggsave(paste0(path2figs,"diversity_stability_relationship_relnonlin.png"), width = 4, height = 4, units = "in", dpi = 120)

ggplot(subset(save_multispp, spprich>0), aes(x=rsd, y=spprich))+
  geom_point()+
  ylab("S")+
  xlab(expression(sigma[R]))+
  theme_classic()+
  theme(axis.title.y = element_text(angle=0, face="italic", size=22),
        axis.title.x = element_text(face="italic", size=22))
ggsave(paste0(path2figs,"diversity_envvar_relationship_relnonlin.png"), width = 3, height = 1.5, units = "in", dpi = 72)


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
  # stat_summary(fun.y = "mean", fill = "coral", size = 4, geom = "point", shape=21)+
  geom_point(shape=21, color="black", fill="dodgerblue", size=3)+
  stat_smooth(method="lm", color="lightblue", fill="grey", se=FALSE)+
  xlab("Number of Species")+
  ylab("Variability of Total Community Biomass (CV)")+
  theme_few()
ggsave(paste0(path2figs,"diversity_stability_relationship_storage_coexist.png"), width = 4, height = 4, units = "in", dpi = 120)

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

ggplot(save_multispp, aes(x=spprich, y=cv))+
  geom_point(shape=21, color="black", fill="coral", size=3)+
  # stat_summary(fun.y = "mean", fill = "coral", size = 5, geom = "point", shape=21)+
  stat_smooth(method="lm", color="orange", fill="grey", se=FALSE)+
  scale_x_continuous(limits=c(1,4))+
  xlab("Number of Species")+
  ylab("Variability of Total Community Biomass (CV)")+
  theme_few()
ggsave(paste0(path2figs,"diversity_stability_relationship_relnonlin_coexist.png"), width = 4, height = 4, units = "in", dpi = 120)

