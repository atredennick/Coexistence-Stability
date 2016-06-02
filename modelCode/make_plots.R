##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

###  PLOT SIMULATION RESULTS ###

# Clear the memory and close the plot windows
rm(list=ls())
for(i in dev.list()) dev.off()

# Check that we are in the correct directory
testdir <- system("ls storageeffect_sims.R", ignore.stderr=TRUE)
if(testdir != 0) {
  cat("\n\nWARNING: Make sure that you are in the correct directory before running the script!\n\n\n")
}

####
####  LOAD LIBRARIES
####
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plyr)
library(reshape2)
library(synchrony)

####
#### INITIALIZATIONS 
####
# Select path to the results and figures
path2results <- "../simulationResults/"
path2figs <- "../manuscript/components/"

seasons_to_exclude <- 500

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
####  PLOT STORAGE EFFECT RESULTS -- CONSTANT RESOURCE
####
## Define vectors of parameters to vary
n_sig_e <- 11 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
n_rho <- 11 # Number of seasonal standard deviation levels
rho_vec <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vector
##  Create matrix with all possible combinations of varying parameters
storage_effect_varvars <- expand.grid(sig_e_vec, rho_vec )
names(storage_effect_varvars) <- c("sigE", "rho")
storage_effect_varvars <- unique(storage_effect_varvars)

storage_effect_sims <- readRDS(paste0(path2results,"storage_effect_equilibrium_runs.RDS"))
save_seasons <- list()
for(i in 1:length(storage_effect_sims)){
  tmp <- as.data.frame(storage_effect_sims[[i]])
  names(tmp) <- c("D1", "D2", "N1", "N2", "R")
  tmp$rho <- storage_effect_varvars[i,"rho"]
  tmp$sigE <- storage_effect_varvars[i,"sigE"]
  tmp$simnum <- i
  tmp$timestep <- 1:nrow(tmp)
  tmp_out <- tmp[seasons_to_exclude:nrow(tmp), ]
  save_seasons <- rbind(save_seasons, tmp_out)
}

# full_melt <- melt(save_seasons, id.vars = c("rho", "sigE", "simnum", "timestep"))
# ggplot(subset(full_melt, variable %in% c("N1", "N2")))+
#   geom_line(aes(x=timestep, y=value, color=variable))+
#   facet_grid(rho~sigE)

##  Calculate CV of Total Community Biomass
save_seasons$total_biomass <- with(save_seasons, N1+N2)
se_community_biomass_cv <- ddply(save_seasons, .(rho, sigE, simnum), summarise,
                                 cv_biomass = sd(total_biomass)/mean(total_biomass))


##  Calculate Invasion Growth Rate
storage_effect_invasions <- readRDS(paste0(path2results,"storage_effect_invasion_runs.RDS"))
se_invasion_growth_rate <- list()
for(i in 1:length(storage_effect_sims)){
  tmp <- storage_effect_invasions[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),2] / 1)
  meantmpr <- mean(tmpr)
  tmpdf <- data.frame(rho=storage_effect_varvars[i,"rho"], sigE=storage_effect_varvars[i,"sigE"], growth_rate=meantmpr)
  se_invasion_growth_rate <- rbind(se_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
se_plot_dat <- merge(se_invasion_growth_rate, se_community_biomass_cv)
se_plot_dat <- subset(se_plot_dat, sigE>0)

##  Make plots
cc_se <- get_colors(length(unique(se_plot_dat$rho)))
se_cv_plot <- ggplot(se_plot_dat, aes(x=sigE, y=cv_biomass, color=as.factor(rho)))+
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

se_gr_plot <- ggplot(se_plot_dat, aes(x=sigE, y=growth_rate, color=as.factor(rho)))+
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

png(paste0(path2figs,"storage_effect_factorial_plots.png"),
    width = 3.5, height=7, units = "in", res=72)
grid_arrange_shared_legend(se_cv_plot, se_gr_plot)
dev.off()


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
for(i in 1:length(storage_effect_sims)){
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
for(i in 1:length(storage_effect_sims)){
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
both_varvars <- matrix(c(-0.5,1,2,2.5,5,5,20,
                -0.5,5,5,20,5,5,20,
                0.5,1,2,2.5,5,5,20,
                0.5,5,5,20,5,5,20), ncol=7, byrow = TRUE)
colnames(both_varvars) <- c("rho", "r1", "a1", "b1", "r2", "a2", "b2")

##  Plot the resource uptake curves
mycols <- c("red", "blue")

mytheme <- theme(axis.text=element_text(size=6),
                 axis.title=element_text(size=8),
                 legend.text=element_text(size=6),
                 legend.key.size=unit(0.3, "cm"))
parms <- list(
  r = c(1,5),         # max growth rate for each species
  alpha = c(2,5),     # rate parameter for Hill function 
  beta = c(2.5,20)    # shape parameter for Hill function
)
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
rel_on_plot <- ggplot()+
  geom_line(data=uptake, aes(x=resource, y=uptake, color=as.factor(species), linetype=as.factor(species)))+
  theme_few()+
  scale_color_manual(values=mycols, name="", labels=c("Species A", "Species B"))+
  guides(size=FALSE, linetype=FALSE)+
  scale_y_continuous(limits=c(0,5))+
  xlab("")+
  ylab("")+
  guides(color=FALSE)+
  mytheme
ggsave("../manuscript/components/Ruptake_relnonlin.png", plot = rel_on_plot, 
       width = 2, height = 1, units = "in", dpi = 72)

# Start low
parms <- list(
  r = c(1,1),         # max growth rate for each species
  alpha = c(2,2),     # rate parameter for Hill function 
  beta = c(2.5,2.5)    # shape parameter for Hill function
)
R <- seq(0,100,1)
out_r <- matrix(ncol=2, nrow=length(R))
alpha <- parms$alpha
beta <- parms$beta
for(i in 1:nrow(out_r)){
  out_r[i,1] <- uptake_R(parms$r[1], R[i], alpha[1], beta[1])
  out_r[i,2] <- uptake_R(parms$r[2], R[i], alpha[2], beta[2])
}
uptake_lo <- data.frame(species=rep(c(1:2), each=nrow(out_r)),
                     resource=rep(R, times=2),
                     uptake=c(out_r[,1], out_r[,2]))
start_lo_plot <- ggplot()+
  geom_line(data=uptake_lo, aes(x=resource, y=uptake, color=as.factor(species), linetype=as.factor(species)))+
  theme_few()+
  scale_color_manual(values=mycols, name="", labels=c("Species A", "Species B"))+
  guides(size=FALSE, linetype=FALSE)+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits=c(0,5))+
  guides(color=FALSE)+
  mytheme
ggsave("../manuscript/components/start_lo_plot.png", plot = start_lo_plot, 
       width = 2, height = 1, units = "in", dpi = 72)


# Start med
parms <- list(
  r = c(3,3),         # max growth rate for each species
  alpha = c(5,5),     # rate parameter for Hill function 
  beta = c(20,20)    # shape parameter for Hill function
)
R <- seq(0,100,1)
out_r <- matrix(ncol=2, nrow=length(R))
alpha <- parms$alpha
beta <- parms$beta
for(i in 1:nrow(out_r)){
  out_r[i,1] <- uptake_R(parms$r[1], R[i], alpha[1], beta[1])
  out_r[i,2] <- uptake_R(parms$r[2], R[i], alpha[2], beta[2])
}
uptake_med <- data.frame(species=rep(c(1:2), each=nrow(out_r)),
                        resource=rep(R, times=2),
                        uptake=c(out_r[,1], out_r[,2]))
start_med_plot <- ggplot()+
  geom_line(data=uptake_med, aes(x=resource, y=uptake, color=as.factor(species), linetype=as.factor(species)))+
  theme_few()+
  scale_color_manual(values=mycols, name="", labels=c("Species A", "Species B"))+
  guides(size=FALSE, linetype=FALSE)+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits=c(0,5))+
  guides(color=FALSE)+
  mytheme
ggsave("../manuscript/components/start_med_plot.png", plot = start_med_plot, 
       width = 2, height = 1, units = "in", dpi = 72)


# Start high
parms <- list(
  r = c(5,5),         # max growth rate for each species
  alpha = c(5,5),     # rate parameter for Hill function 
  beta = c(20,20)    # shape parameter for Hill function
)
R <- seq(0,100,1)
out_r <- matrix(ncol=2, nrow=length(R))
alpha <- parms$alpha
beta <- parms$beta
for(i in 1:nrow(out_r)){
  out_r[i,1] <- uptake_R(parms$r[1], R[i], alpha[1], beta[1])
  out_r[i,2] <- uptake_R(parms$r[2], R[i], alpha[2], beta[2])
}
uptake_hi <- data.frame(species=rep(c(1:2), each=nrow(out_r)),
                         resource=rep(R, times=2),
                         uptake=c(out_r[,1], out_r[,2]))
start_hi_plot <- ggplot()+
  geom_line(data=uptake_hi, aes(x=resource, y=uptake, color=as.factor(species), linetype=as.factor(species)))+
  theme_few()+
  scale_color_manual(values=mycols, name="", labels=c("Species A", "Species B"))+
  guides(size=FALSE, linetype=FALSE)+
  xlab("")+
  ylab("")+
  scale_y_continuous(limits=c(0,5))+
  guides(color=FALSE)+
  mytheme
ggsave("../manuscript/components/start_hi_plot.png", plot = start_hi_plot, 
       width = 2, height = 1, units = "in", dpi = 72)


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
  save_seasons_both <- rbind(save_seasons_both, tmp_out)
}

##  Calculate CV of Total Community Biomass
save_seasons_both$total_biomass <- with(save_seasons_both, N1+N2)
both_community_biomass_cv <- ddply(save_seasons_both, .(rho, rel_nonlin), summarise,
                                 cv_biomass = sd(total_biomass)/mean(total_biomass))

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
  both_invasion_growth_rate <- rbind(both_invasion_growth_rate, tmpdf)
}

both_together <- merge(both_community_biomass_cv, both_invasion_growth_rate)

both1 <- ggplot(both_together, aes(x=as.factor(rho),y=cv_biomass, group=rel_nonlin))+
  geom_line()+
  geom_point(size=4,color="white")+
  geom_point(size=4, aes(shape=rel_nonlin))+
  ylab("CV of community biomass")+
  xlab(expression(rho))+
  theme_few()+
  scale_shape_manual(values=c(1,19), name="", labels=c("RNL absent", "RNL present"))+
  theme(legend.position=c(0.75,0.15))+
  theme(legend.text=element_text(size=10),
        legend.key.size=unit(0.5, "cm"))+
  guides(shape = guide_legend(override.aes = list(size=2)))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.3, "cm"))

both2 <- ggplot(both_together, aes(x=as.factor(rho),y=growth_rate, group=rel_nonlin))+
  geom_hline(aes(yintercept=0))+
  geom_line()+
  geom_point(size=4,color="white")+
  geom_point(size=4, aes(shape=rel_nonlin))+
  ylab("Log invasion growth rate")+
  xlab(expression(rho))+  theme_few()+
  scale_shape_manual(values=c(1,19), name="", labels=c("Relative nonlinearity absen", "Relative nonlinearity absent"))+
  guides(shape=FALSE)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15),
        legend.text=element_text(size=6),
        legend.key.size=unit(0.3, "cm"))

png(paste0(path2figs,"storage_effect_relative_nonlinearity_anova_plot.png"),
    width = 3.5, height=7, units = "in", res=72)
grid.arrange(both1, both2)
dev.off()
