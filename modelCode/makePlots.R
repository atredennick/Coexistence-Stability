##============================================================================##
##                                                                            ##
## A. Tredennick, P. Adler, and F. Adler                                      ##             
## "How Fluctuation-Dependent Coexistence Mechanisms Affect the Temporal      ##
##  Stability of Ecosystem Function"                                          ##
##                                                                       2016 ##
##============================================================================##

##  PLOT SIMULATION RESULTS

# Clear the memory and close the plot windows
rm(list=ls())

library(ggplot2)
library(ggthemes)
library(gridExtra)
library(plyr)
library(reshape2)
library(synchrony)

####
#### INITIALIZATIONS 
####
# for(i in dev.list()) dev.off()

# Check that we are in the correct directory
testdir <- system("ls storageeffect_sims.R", ignore.stderr=TRUE)
if(testdir!= 0){
  cat("\n\nWARNING: Make sure that you are in the correct directory before running the script!\n\n\n")
}

# Select path to the results and figures
path2results <- "../simulationResults/"
path2figs <- "../manuscript/components/"

# Recreate simulation grid
nrho <- 11
rholist <- pretty(seq(-1, 1, length.out=nrho), nrho)
nsd <- 11
rsdlist <- pretty(seq(0, 1, length.out=nsd), nsd)
storage_effect_varvars <- expand.grid(rholist,rsdlist)
names(storage_effect_varvars) <- c("rho", "Rsd")



####
####  PLOT STORAGE EFFECT RESULTS
####

##  1. Species asynchrony vs. environmental cue correlation
sims <- readRDS(paste0(path2results,"storage_effect_sims.RDS"))

#take out first couple seasons
seasons_to_exclude <- 20
save_seasons <- data.frame(time=NA, D1=NA, D2=NA, N1=NA, N2=NA, R=NA, Rstart=NA,
                           season=NA, rho=NA, rsd=NA, simnum=NA)
for(i in 1:length(sims)){
  tmp <- sims[[i]]
  tmp <- tmp[2:nrow(tmp),]
  tmp$rho <- storage_effect_varvars[i,"rho"]
  tmp$rsd <- storage_effect_varvars[i,"Rsd"]
  tmp$simnum <- i
  save_seasons <- rbind(save_seasons, tmp)
}
save_seasons <- save_seasons[2:nrow(save_seasons),]

# Analyze average biomass over the season
mean_of_season <- ddply(save_seasons, .(season, rho, rsd, simnum), summarise,
                        mean_N1 = mean(N1),
                        mean_N2 = mean(N2),
                        mean_D1 = mean(D1),
                        mean_D2 = mean(D2))

tmp.sync <- numeric(length(sims))
for(i in 1:length(sims)){
  tmp <- subset(mean_of_season, simnum==i)
  tmp.sync[i] <- as.numeric(community.sync(tmp[,c("mean_N1", "mean_N2")])[1])
}

syncdf <- data.frame(rho=storage_effect_varvars$rho,
                     rsd=storage_effect_varvars$Rsd,
                     async=1-tmp.sync)

lowcol <- "black"
highcol <- "grey"

sync.strg <- ggplot()+
  geom_line(data=syncdf, aes(x=rho, y=async, color=rsd, group=as.factor(rsd)))+
  geom_point(data=syncdf, aes(x=rho, y=async, color=rsd, group=as.factor(rsd)),size=2)+
  theme_few()+
  xlab(bquote(rho))+
  ylab("Species Asynchrony")+
  # scale_color_continuous(name=bquote(sigma[R]))+
  scale_y_continuous(limits=c(0,1))+
  scale_colour_gradient(low=lowcol, high=highcol, name=bquote(sigma[R]))+
  geom_text(aes(x=1, y=1, label="a"))

seasonal_total <- ddply(mean_of_season, .(season, rho, rsd, simnum), summarise,
                        total_biomass=sum(mean_N1, mean_N2))
stability <- ddply(seasonal_total, .(rho, rsd, simnum), summarise,
                   stable = sd(total_biomass)/mean(total_biomass),
                   mubiom = mean(total_biomass),
                   sdbiom = sd(total_biomass))

# Get correlation of stable coexistenc (rho) and ecosystem stability (CV)
cor((1-stability$rho), stability$stable, method = "spearman")

stab.strg <- ggplot()+
  geom_line(data=stability, aes(x=rho, y=stable, color=rsd, group=as.factor(rsd)))+
  geom_point(data=stability, aes(x=rho, y=stable, color=rsd, group=as.factor(rsd)), size=2)+
  theme_few()+
  xlab(bquote(rho))+
  ylab("CV of Total Biomass")+
  scale_colour_gradient(low=lowcol, high=highcol, name=bquote(sigma[R]))+
  scale_y_continuous(limits=c(0,1))+
  geom_text(aes(x=1, y=1, label="b"))

mu.strg <- ggplot()+
  geom_line(data=stability, aes(x=rho, y=mubiom, color=rsd, group=as.factor(rsd)))+
  geom_point(data=stability, aes(x=rho, y=mubiom, color=rsd, group=as.factor(rsd)), size=2)+
  theme_few()+
  xlab(bquote(rho))+
  ylab("Mean Total Biomass")+
  scale_colour_gradient(low=lowcol, high=highcol, name=bquote(sigma[R]))+
  scale_y_continuous(limits=c(0,100))+
  geom_text(aes(x=1, y=100, label="c"))

sd.strg <- ggplot()+
  geom_line(data=stability, aes(x=rho, y=sdbiom, color=rsd, group=as.factor(rsd)))+
  geom_point(data=stability, aes(x=rho, y=sdbiom, color=rsd, group=as.factor(rsd)), size=2)+
  theme_few()+
  xlab(bquote(rho))+
  ylab("SD Total Biomass")+
  scale_colour_gradient(low=lowcol, high=highcol, name=bquote(sigma[R]))+
  scale_y_continuous(limits=c(0,100))+
  geom_text(aes(x=1, y=100, label="d"))

# strg_results_fig <- grid.arrange(sync.strg, stab.strg, ncol=2)
png(paste0(path2figs,"storage_effect_synchstability.png"), width = 8.5, height = 4, 
    units="in", res=100)
strg_results_fig <- grid.arrange(sync.strg, stab.strg, mu.strg, sd.strg, ncol=2, nrow=2)
print(strg_results_fig)
dev.off()

strg_stability <- stability

####
####  RELATIVE NONLINEARITY PLOTS
####
# Recreate simulation grid
nrho <- 21
rholist <- rep(1, nrho)
nsd <- 21
rsdlist <- pretty(seq(0, 1, length.out=nsd), nsd)
relnonlin_varvars <- expand.grid(rholist,rsdlist)
names(relnonlin_varvars) <- c("rho", "Rsd")
relnonlin_varvars <- unique(relnonlin_varvars)

##  1. Species asynchrony vs. environmental cue correlation
sims <- readRDS(paste0(path2results,"relative_nonlinearity_sims.RDS"))

#take out first couple seasons
seasons_to_exclude <- 20
save_seasons <- data.frame(time=NA, D1=NA, D2=NA, N1=NA, N2=NA, R=NA, Rstart=NA,
                           season=NA, rho=NA, rsd=NA, simnum=NA)
for(i in 1:length(sims)){
  tmp <- sims[[i]]
  tmp <- tmp[2:nrow(tmp),]
  tmp$rho <- relnonlin_varvars[i,"rho"]
  tmp$rsd <- relnonlin_varvars[i,"Rsd"]
  tmp$simnum <- i
  save_seasons <- rbind(save_seasons, tmp)
}
save_seasons <- save_seasons[2:nrow(save_seasons),]

# Analyze average biomass over the season
mean_of_season <- ddply(save_seasons, .(season, rho, rsd, simnum), summarise,
                        mean_N1 = mean(N1),
                        mean_N2 = mean(N2),
                        mean_D1 = mean(D1),
                        mean_D2 = mean(D2))

tmp.sync <- numeric(length(sims))
for(i in 1:length(sims)){
  tmp <- subset(mean_of_season, simnum==i)
  tmp.sync[i] <- as.numeric(community.sync(tmp[,c("mean_N1", "mean_N2")])[1])
}

syncdf <- data.frame(rho=relnonlin_varvars$rho,
                     rsd=relnonlin_varvars$Rsd,
                     async=1-tmp.sync)

sync.strg <- ggplot()+
  geom_line(data=syncdf, aes(x=rsd, y=async))+
  geom_point(data=syncdf, aes(x=rsd, y=async),size=2)+
  theme_few()+
  xlab(bquote(sigma[R]))+
  ylab("Species Asynchrony")+
  scale_y_continuous(limits=c(0,1))+
  geom_text(aes(x=1, y=1, label="a"))

seasonal_total <- ddply(mean_of_season, .(season, rho, rsd, simnum), summarise,
                        total_biomass=sum(mean_N1, mean_N2))
stability <- ddply(seasonal_total, .(rho, rsd, simnum), summarise,
                   stable = sd(total_biomass)/mean(total_biomass),
                   mubiom = mean(total_biomass),
                   sdbiom = sd(total_biomass))

# Get correlation of stable coexistenc (rho) and ecosystem stability (CV)
cor((stability$rsd), stability$stable, method = "spearman")

stab.strg <- ggplot()+
  geom_line(data=stability, aes(x=rsd, y=stable))+
  geom_point(data=stability, aes(x=rsd, y=stable), size=2)+
  theme_few()+
  xlab(bquote(sigma[R]))+
  ylab("CV of Total Biomass")+
  scale_y_continuous(limits=c(0,1))+
  geom_text(aes(x=1, y=1, label="b"))

mu.strg <- ggplot()+
  geom_line(data=stability, aes(x=rsd, y=mubiom))+
  geom_point(data=stability, aes(x=rsd, y=mubiom), size=2)+
  theme_few()+
  xlab(bquote(sigma[R]))+
  ylab("Mean Total Biomass")+
  scale_y_continuous(limits=c(0,100))+
  geom_text(aes(x=1, y=100, label="c"))

sd.strg <- ggplot()+
  geom_line(data=stability, aes(x=rsd, y=sdbiom))+
  geom_point(data=stability, aes(x=rsd, y=sdbiom), size=2)+
  theme_few()+
  xlab(bquote(sigma[R]))+
  ylab("SD Total Biomass")+
  scale_y_continuous(limits=c(0,100))+
  geom_text(aes(x=1, y=100, label="d"))

# strg_results_fig <- grid.arrange(sync.strg, stab.strg, ncol=2)
png(paste0(path2figs,"relative_nonlinearity_synchstability.png"), width = 8.5, height = 4, 
    units="in", res=100)
strg_results_fig <- grid.arrange(sync.strg, stab.strg, mu.strg, sd.strg, ncol=2, nrow=2)
print(strg_results_fig)
dev.off()

relnonlin_stability <- stability



####
####  PLOT INVASION GROWTH RATE V. CV OF BIOMASS
####
##  Relative Nonlinearity
# Recreate simulation grid
nrho <- 21
rholist <- rep(1, nrho)
nsd <- 21
rsdlist <- pretty(seq(0, 1, length.out=nsd), nsd)
relnonlin_varvars <- expand.grid(rholist,rsdlist)
names(relnonlin_varvars) <- c("rho", "Rsd")
relnonlin_varvars <- unique(relnonlin_varvars)

sims <- readRDS(paste0(path2results,"relative_nonlinearity_invasion_sims.RDS"))

#take out first couple seasons
seasons_to_exclude <- 20
save_seasons <- data.frame(time=NA, D1=NA, D2=NA, N1=NA, N2=NA, R=NA, Rstart=NA,
                           season=NA, rho=NA, rsd=NA, simnum=NA)
for(i in 1:length(sims)){
  tmp <- sims[[i]]
  tmp <- tmp[2:nrow(tmp),]
  tmp$rho <- relnonlin_varvars[i,"rho"]
  tmp$rsd <- relnonlin_varvars[i,"Rsd"]
  tmp$simnum <- i
  save_seasons <- rbind(save_seasons, tmp)
}
save_seasons <- save_seasons[2:nrow(save_seasons),]

# Analyze average biomass over the season
mean_of_season <- ddply(save_seasons, .(season, rho, rsd, simnum), summarise,
                        mean_N1 = mean(N1),
                        mean_N2 = mean(N2),
                        mean_D1 = mean(D1),
                        mean_D2 = mean(D2))
invasion_growth_rates <- mean_of_season[,c("season", "rho", "rsd", "mean_N1")]
invasion_growth_rates$igr <- log(invasion_growth_rates$mean_N1) - log(0.01)
mean_invasion_growth <- ddply(invasion_growth_rates, .(rho, rsd), summarise,
                              avg_igr = mean(igr))
igr_df_rel <- data.frame(strength_coexist = mean_invasion_growth$avg_igr, 
                     cv_biomass = relnonlin_stability$stable)
cor(igr_df_rel$strength_coexist, igr_df_rel$cv_biomass, method = "spearman")

##  Storage Effect
sims <- readRDS(paste0(path2results,"storage_effect_invasion_sims.RDS"))

#take out first couple seasons
seasons_to_exclude <- 20
save_seasons <- data.frame(time=NA, D1=NA, D2=NA, N1=NA, N2=NA, R=NA, Rstart=NA,
                           season=NA, rho=NA, rsd=NA, simnum=NA)
for(i in 1:length(sims)){
  tmp <- sims[[i]]
  tmp <- tmp[2:nrow(tmp),]
  tmp$rho <- storage_effect_varvars[i,"rho"]
  tmp$rsd <- storage_effect_varvars[i,"Rsd"]
  tmp$simnum <- i
  save_seasons <- rbind(save_seasons, tmp)
}
save_seasons <- save_seasons[2:nrow(save_seasons),]

# Analyze average biomass over the season
mean_of_season <- ddply(save_seasons, .(season, rho, rsd, simnum), summarise,
                        mean_N1 = mean(N1),
                        mean_N2 = mean(N2),
                        mean_D1 = mean(D1),
                        mean_D2 = mean(D2))
invasion_growth_rates <- mean_of_season[,c("season", "rho", "rsd", "mean_N2")]
invasion_growth_rates$igr <- log(invasion_growth_rates$mean_N2) - log(0.01)
mean_invasion_growth <- ddply(invasion_growth_rates, .(rho, rsd), summarise,
                              avg_igr = mean(igr))
igr_df_strg <- data.frame(strength_coexist = mean_invasion_growth$avg_igr, 
                          cv_biomass = strg_stability$stable,
                          rho = strg_stability$rho,
                          rsd = strg_stability$rsd)
cor(igr_df_strg$strength_coexist, igr_df_strg$cv_biomass, method = "spearman")

igr_df_strg <- igr_df_strg[-55,]
ggplot(igr_df_strg, aes(x=strength_coexist, y=cv_biomass, color=rho))+
  geom_point()

strength_plot <- ggplot(igr_df, aes(x=strength_coexist, y=cv_biomass))+
  geom_point(size=4)+
  xlab("Strength of Coexistence")+
  ylab("CV of Total Community Biomass")+
  theme_few()

png(paste0(path2figs,"strength_v_cv.png"), width = 5, height = 4, 
    units="in", res=100)
print(strength_plot)
dev.off()


####
####  COMPARE CV UNDER EACH COEXISTENCE MECHANISM
####
relnonlin_stability <- data.frame(rho = rep(pretty(seq(-1, 1, length.out=11), 11), each=11),
                                  rsd = rep(relnonlin_stability$rsd, times=11),
                                  simnum = strg_stability$simnum,
                                  stable = rep(relnonlin_stability$stable, times=11))
relnonlin_stability$type <- "relnon"
strg_stability$type <- "storage"
combo_stability <- rbind(relnonlin_stability, strg_stability[,c("rho", "rsd", "simnum", "stable", "type")])

mycolors <- c("#9D6188","#97A861")
compare_plot <- ggplot(data=combo_stability, aes(x=rsd, y=stable, color=type))+
  geom_point()+
  geom_line()+
  facet_wrap("rho")+
  xlab(bquote(sigma[R]))+
  ylab("CV of Total Community Biomass")+
  scale_y_continuous(limits=c(0,1))+
  scale_color_manual(values=mycolors, labels=c("Relative Nonlinearity", "Storage Effect"), name="")+
  theme_few()+
  theme(axis.text.x=element_text(angle=45, hjust=1))

png(paste0(path2figs,"compare_plots.png"), width = 8.5, height = 5, 
    units="in", res=100)
print(compare_plot)
dev.off()


## Plot the difference between coexistence mechanisms
# wide_stability <- dcast(combo_stability, rho+rsd~type, value.var = "stable")
# wide_stability$diff <- with(wide_stability, relnon-storage)
# ggplot(wide_stability, aes(x=rsd, y=diff, color=rho, group=rho))+
#   # geom_line()+
#   geom_point()

