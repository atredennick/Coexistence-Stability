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
if(testdir != 0) {
  cat("\n\nWARNING: Make sure that you are in the correct directory before running the script!\n\n\n")
}

# Select path to the results and figures
path2results <- "../simulationResults/"
path2figs <- "../manuscript/components/"

seasons_to_exclude <- 0


####
####  PLOT STORAGE EFFECT RESULTS
####
## Define vectors of parameters to vary
n_sig_e <- 11 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 5, length.out=n_sig_e), n_sig_e) # Make a pretty vector
n_rho <- 11 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(-1, 1, length.out=n_rho), n_rho) # Make a pretty vectorhttp://127.0.0.1:45503/graphics/plot_zoom_png?width=1200&height=900
##  Create matrix with all possible combinations of varying parameters
storage_effect_varvars <- expand.grid(sig_e_vec, rsd_vec )
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

full_melt <- melt(save_seasons, id.vars = c("rho", "sigE", "simnum", "timestep"))
ggplot(subset(full_melt, variable %in% c("N1", "N2")))+
  geom_line(aes(x=timestep, y=value, color=variable))+
  facet_grid(rho~sigE)

##  Calculate CV of Total Community Biomass
save_seasons$total_biomass <- with(save_seasons, N1+N2)
se_community_biomass_cv <- ddply(save_seasons, .(rho, sigE, simnum), summarise,
                                 cv_biomass = sd(total_biomass)/mean(total_biomass))


##  Calculate Invasion Growth Rate
storage_effect_invasions <- readRDS(paste0(path2results,"storage_effect_invasion_runs.RDS"))
se_invasion_growth_rate <- list()
for(i in 1:length(storage_effect_sims)){
  tmp <- storage_effect_invasions[[i]]
  tmpr <- log(tmp[seasons_to_exclude:nrow(tmp),4] / 1)
  meantmpr <- mean(tmpr)
  tmpdf <- data.frame(rho=storage_effect_varvars[i,"rho"], sigE=storage_effect_varvars[i,"sigE"], growth_rate=meantmpr)
  se_invasion_growth_rate <- rbind(se_invasion_growth_rate, tmpdf)
}

##  Merge Growth Rates and Biomass CV; Plot
se_plot_dat <- merge(se_invasion_growth_rate, se_community_biomass_cv)
se_plot_dat <- subset(se_plot_dat, sigE>0)
ggplot(se_plot_dat, aes(x=growth_rate, y=cv_biomass, color=rho))+
  geom_point()+
  stat_smooth(method="lm", color="black", se=TRUE)+
  facet_wrap("sigE", ncol=10)+
  theme_few()



out <- by(data = se_plot_dat, INDICES = list(se_plot_dat$sigE), FUN = function(x) {
  model <- lm(cv_biomass ~ growth_rate, data = x)
  data.frame(
    intercept = coef(model)["(Intercept)"], 
    slope = coef(model)["growth_rate"], 
    id = unique(x$sigE)
  )
})
all_slopes <- do.call("rbind", out)

ggplot(all_slopes, aes(x=id, y=slope))+
  geom_point(size=5)+
  stat_smooth(se=FALSE, size=1)
