##  collate_div+envvar_sims.R

rm(list=ls(all.names = TRUE))


#===========================#
####### PRELIMINARIES #######
#===========================#

####
####  Libraries
####
library("plyr")
library("reshape2")
library("ggplot2")



####
####  Function for calculating CV from list elements
####
get_cv <- function(x){
  seasons_to_exclude <- 500 # always get rid of the first 500 iterations
  x <- as.data.frame(x)
  names(x) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(x))
  tmp_totbiomass <- rowSums(x[seasons_to_exclude:nrow(x),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(x[seasons_to_exclude:nrow(x),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  tmp_out <- c(tmp_cv, tmp_spprich)
}



####
####  Function for collating raw simulation lists
####
collate_sims <- function(in_file, parameter_matrix, save_file){
  sims_list <- readRDS(in_file)                     # read in the large simulation output
  out_cvs <- lapply(sims_list, get_cv)              # apply the function
  out_cvs_df <- do.call(rbind.data.frame, out_cvs)  # convert to data frame
  colnames(out_cvs_df) <- c("cv", "spprichness")    # give the columns names
  cvs_params <- cbind(parameter_matrix, out_cvs_df) # combine with parameter values
  saveRDS(cvs_params, file = save_file)             # save the file
  # file.remove(in_file)                              # delete the large, raw file
}



####
####  Function for returning the slope and 95% CI
####
get_slope <- function(df_subset) {
  mod <- lm(log(cv) ~ log(sigE), data = df_subset)
  slope <- c(coef(mod)[2], confint(mod)[2,1], confint(mod)[2,2])
  return(slope)
}

get_slope_rnonlin <- function(df_subset) {
  mod <- lm(log(cv) ~ log(Rsd_annual), data = df_subset)
  slope <- c(coef(mod)[2], confint(mod)[2,1], confint(mod)[2,2])
  return(slope)
}



####
####  My plotting theme (for ggplot2)
####
my_theme <- theme(legend.title=element_text(size=8, face="bold"),
                  legend.text=element_text(size=8),
                  legend.background = element_rect(colour = "grey45", size=0.5),
                  legend.key = element_blank(),
                  legend.key.size = unit(0.3, "cm"),
                  panel.grid.major = element_line(colour = "grey25", linetype = "dotted"),
                  panel.grid.minor = element_line(colour = "white", linetype = "dotted"),
                  strip.background = element_blank())



#========================================#
####### STORAGE EFFECT SIMULATIONS #######
#========================================#

####
####  Compare symmetric and asymmetric competition
####

##  Recreate parameter matrix for simulations
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))
n_sig_e <- 15 # Number of cue variance levels
sig_e_vec <- pretty(seq(0.1, 2, length.out=n_sig_e), n_sig_e) # Make a pretty vector
rho <- as.matrix(c(-1,0,1))
prm <- expand.grid(as.matrix(sig_e_vec), rho, 1:4)
DNR_repped <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
colnames(prm) <- c("sigE", "rho", "dnr_id")
colnames(DNR_repped) <- colnames(DNR)
prm <- cbind(prm, DNR_repped)
prm <- subset(prm, select = -c(dnr_id))
alphas <- rbind(c(alpha1=0.5, alpha2=0.495, alpha3=0.49, alpha4=0.485),
                c(alpha1=0.5, alpha2=0.49, alpha3=0.48, alpha4=0.47))
prm_repped <- do.call("rbind", replicate(2, prm,  simplify = FALSE))
alphas_repped <- matrix(rep(alphas, each=nrow(prm)), ncol=ncol(alphas))
colnames(alphas_repped) <- colnames(alphas)
parameter_matrix <- cbind(prm_repped, alphas_repped)


##  Collate the results, if necessary, otherwise read in collated results
strg_comp_large_file <- "../simulationResults/storageeffect_div+envvar_stability_varycomp.RDS"
strg_comp_file <- "../simulationResults/storageeffect_div+envvar_cv_varycomp_collated.RDS"
if(file.exists(strg_comp_file) == FALSE) {
  collate_sims(in_file = strg_comp_large_file, 
               parameter_matrix = parameter_matrix, 
               save_file = strg_comp_file)
}


##  Read in collated results and do some housekeeping columns
strg_comp_cvs <- readRDS(strg_comp_file)
strg_comp_cvs$comp <- "1symmetric"
strg_comp_cvs[which(strg_comp_cvs$alpha4==0.47), "comp"] <- "2asymmetric"

##  Summarise by taking mean CV at different realized richness
strg_comp_cvs_mean <- ddply(strg_comp_cvs, .(sigE, rho, spprichness, comp), 
                            summarise,
                            cv = mean(cv))

## Find species first occurences for lines
firstones <- ddply(strg_comp_cvs_mean, 
                   .(rho, comp, spprichness), 
                   summarise, 
                   sigE=min(sigE),
                   cv=cv[which.min(sigE)])

##  Plot and save
ggplot()+
  geom_segment(data=firstones, aes(x=sigE, y=cv, xend=sigE, yend=0, color=as.factor(spprichness)), size=0.5, alpha=0.7)+
  geom_line(data=strg_comp_cvs_mean, aes(x=sigE, y=cv, color=as.factor(spprichness)), size=0.7, alpha=0.7)+
  geom_point(data=strg_comp_cvs_mean, aes(x=sigE, y=cv, color=as.factor(spprichness)), size=2, alpha=1)+
  geom_point(data=strg_comp_cvs_mean, aes(x=sigE, y=cv, color=as.factor(spprichness)), size=2, alpha=1, shape=1, color="grey15")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(comp~rho)+
  theme_bw()+
  my_theme+
  theme(legend.position = c(0.94, 0.155))+
  theme(legend.title=element_text(size=6, face="bold"),
        legend.text=element_text(size=6),
        legend.background = element_rect(colour = "grey45", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, "cm"))+
  guides(colour = guide_legend(override.aes = list(size=1)))
ggsave(filename = "../manuscript/components/storage_effect_div+envar_varycomp.png", width = 133, height=80, units = "mm", dpi = 600)


## Log-log plot and save
ggplot(strg_comp_cvs_mean, aes(x=log(sigE), y=log(cv), color=as.factor(spprichness)))+
  stat_smooth(method="lm", size=0.7, alpha=0.7, se=FALSE)+
  geom_point(size=2, alpha=1)+
  geom_point(size=2, alpha=1, shape=1, color="grey15")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(comp~rho)+
  theme_bw()+
  my_theme+
  theme(legend.position = c(0.95, 0.15))
ggsave(filename = "../manuscript/components/storage_effect_div+envar_varycomp_loglog.png", width = 8, height=4, units = "in", dpi = 82)


##  Get slopes and plot
strg_comp_cvs_slopes <- list()
for(do_comp in unique(strg_comp_cvs_mean$comp)){
  for(do_rho in c(-1,0)){
    tmp1 <- subset(strg_comp_cvs_mean, comp==do_comp & rho==do_rho)
    for(do_spprich in unique(tmp1$spprichness)){
      tmp_df <- subset(tmp1, spprichness==do_spprich)
      tmp_slope <- get_slope(tmp_df)
      tmpout <- data.frame(comp=do_comp, spprich=do_spprich, rho=do_rho, 
                           slope=tmp_slope[1], lowslope=tmp_slope[2], hislope=tmp_slope[3])
      strg_comp_cvs_slopes <- rbind(strg_comp_cvs_slopes, tmpout)
    }
  }
}

ylabel <- bquote("Slope (CV vs." ~ sigma[E] ~ ")")
ggplot( strg_comp_cvs_slopes, aes(x=spprich, y=slope, color=as.factor(comp)))+
  geom_ribbon(aes(ymin=lowslope, ymax=hislope, fill=as.factor(comp)), alpha=0.2, color=NA)+
  geom_point()+
  geom_line()+
  facet_wrap("rho")+
  ylab(ylabel)+
  xlab("Realized species richness")+
  scale_fill_manual(values=c("red", "blue"), name="Competition", labels=c("Symmetric comp.", "Asymmetric comp."))+
  scale_color_manual(values=c("red", "blue"), name="Competition", labels=c("Symmetric comp.", "Asymmetric comp."))+
  theme_bw()+
  my_theme
ggsave("../manuscript/components/storage_effect_div+envar_varycomp_loglog_slopes.png", width = 5, height=2, units = "in", dpi = 120)



####
####  Compare high and low dormant mortality
####
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))
n_sig_e <- 15 # Number of cue variance levels
sig_e_vec <- pretty(seq(0.1, 2, length.out=n_sig_e), n_sig_e) # Make a pretty vector
rho <- as.matrix(c(-1,0,1))
prm <- expand.grid(as.matrix(sig_e_vec), rho, 1:4)
DNR_repped <- matrix(rep(DNR,each=(nrow(prm)/nrow(DNR))),ncol=ncol(DNR))
colnames(prm) <- c("sigE", "rho", "dnr_id")
colnames(DNR_repped) <- colnames(DNR)
prm <- cbind(prm, DNR_repped)
prm <- subset(prm, select = -c(dnr_id))
loweta <- 0.2
hieta <- 0.8
etas <- rbind(c(eta1=loweta, eta2=loweta, eta3=loweta, eta4=loweta),
              c(eta1=hieta, eta2=hieta, eta3=hieta, eta4=hieta))
prm_repped <- do.call("rbind", replicate(2, prm,  simplify = FALSE))
etas_repped <- matrix(rep(etas, each=nrow(prm)), ncol=ncol(etas))
colnames(etas_repped) <- colnames(etas)
parameter_matrix <- cbind(prm_repped, etas_repped)


##  Collate the results, if necessary, otherwise read in collated results
strg_eta_large_file <- "../simulationResults/storageeffect_div+envvar_stability_varyeta.RDS"
strg_eta_file <- "../simulationResults/storageeffect_div+envvar_cv_varyeta_collated.RDS"
if(file.exists(strg_eta_file) == FALSE) {
  collate_sims(in_file = strg_eta_large_file, 
               parameter_matrix = parameter_matrix, 
               save_file = strg_eta_file)
}


##  Read in collated results and do some housekeeping columns
strg_eta_cvs <- readRDS(strg_eta_file)

##  Summarise by taking mean CV at different realized richness
strg_eta_cvs_mean <- ddply(strg_eta_cvs, .(sigE, rho, spprichness, eta1), 
                            summarise,
                            cv = mean(cv))

##  Plot and save
ggplot(strg_eta_cvs_mean, aes(x=sigE, y=cv, color=as.factor(spprichness)))+
  geom_line(size=0.7, alpha=0.7)+
  geom_point(size=2, alpha=1)+
  geom_point(size=2, alpha=1, shape=1, color="grey15")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(eta1~rho, scales="free_y")+
  theme_bw()+
  my_theme+
  theme(legend.position = c(0.95, 0.15))
ggsave(filename = "../manuscript/components/storage_effect_div+envar_varyeta.png", width = 8, height=4, units = "in", dpi = 82)

## Log-log plot and save
ggplot(strg_eta_cvs_mean, aes(x=log(sigE), y=log(cv), color=as.factor(spprichness)))+
  stat_smooth(method="lm", size=0.7, alpha=0.7, se=FALSE)+
  geom_point(size=2, alpha=1)+
  geom_point(size=2, alpha=1, shape=1, color="grey15")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(eta1~rho, scales="free_y")+
  theme_bw()+
  my_theme+
  theme(legend.position = c(0.95, 0.15))
ggsave(filename = "../manuscript/components/storage_effect_div+envar_varyeta_loglog.png", width = 8, height=4, units = "in", dpi = 82)

##  Get slopes and plot
strg_eta_cvs_slopes <- list()
for(do_eta in unique(strg_eta_cvs_mean$eta1)){
  for(do_rho in c(-1,0)){
    for(do_spprich in unique(strg_eta_cvs_mean$spprichness)){
      tmp_df <- subset(strg_eta_cvs_mean, eta1==do_eta & spprichness==do_spprich & rho==do_rho)
      tmp_slope <- get_slope(tmp_df)
      tmpout <- data.frame(eta=do_eta, spprich=do_spprich, rho=do_rho, 
                           slope=tmp_slope[1], lowslope=tmp_slope[2], hislope=tmp_slope[3])
      strg_eta_cvs_slopes <- rbind(strg_eta_cvs_slopes, tmpout)
    }
  }
}

ylabel <- bquote("Slope (CV vs." ~ sigma[E]^2 ~ ")")
ggplot( strg_eta_cvs_slopes, aes(x=spprich, y=slope, color=as.factor(eta)))+
  geom_ribbon(aes(ymin=lowslope, ymax=hislope, fill=as.factor(eta)), alpha=0.2, color=NA)+
  geom_point()+
  geom_line()+
  facet_wrap("rho")+
  ylab(ylabel)+
  xlab("Realized species richness")+
  scale_fill_manual(values=c("red", "blue"), name=expression(eta), labels=c("Symmetric comp.", "Asymmetric comp."))+
  scale_color_manual(values=c("red", "blue"), name=expression(eta), labels=c("Symmetric comp.", "Asymmetric comp."))+
  theme_bw()+
  my_theme
ggsave("../manuscript/components/storage_effect_div+envar_varyeta_loglog_slopes.png", width = 5, height=2, units = "in", dpi = 120)








#================================================#
####### RELATIVE NONLINEARITY SIMULATIONS #######
#===============================================#

####
####  Stable to unstable and unstable to stable species additions
####
##  Recreate parameter matrix -- stable to unstable
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))

## Define vectors of parameters to vary
n_rsd <- 25 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0.1, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
prm <- as.data.frame(rsd_vec)
colnames(prm) <- "Rsd_annual"
# Add in variable parameters to form parameter matrix
DNR_repped <- do.call("rbind", replicate(nrow(prm), DNR,  simplify = FALSE))
prm_repped <- do.call("rbind", replicate(nrow(DNR), prm,  simplify = FALSE))
parameter_matrix <- cbind(DNR_repped, prm_repped)

##  Collate the results, if necessary, otherwise read in collated results
rnonlin_stable2unstable_large_file <- "../simulationResults/relnonlin_div+envvar_cv_stable2unstable.RDS"
rnonlin_stable2unstable_file <- "../simulationResults/relnonlin_div+envvar_cv_stable2unstable_collated.RDS"
if(file.exists(rnonlin_stable2unstable_file) == FALSE) {
  collate_sims(in_file = rnonlin_stable2unstable_large_file, 
               parameter_matrix = parameter_matrix, 
               save_file = rnonlin_stable2unstable_file)
}


##  Recreate parameter matrix -- unstable to stable
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(0,1,1,1),N=c(0,1,1,1),R=20),
             c(D=c(0,0,1,1),N=c(0,0,1,1),R=20),
             c(D=c(0,0,0,1),N=c(0,0,0,1),R=20))

## Define vectors of parameters to vary
n_rsd <- 25 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0.1, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
prm <- as.data.frame(rsd_vec)
colnames(prm) <- "Rsd_annual"
# Add in variable parameters to form parameter matrix
DNR_repped <- do.call("rbind", replicate(nrow(prm), DNR,  simplify = FALSE))
prm_repped <- do.call("rbind", replicate(nrow(DNR), prm,  simplify = FALSE))
parameter_matrix <- cbind(DNR_repped, prm_repped)

##  Collate the results, if necessary, otherwise read in collated results
rnonlin_unstable2stable_large_file <- "../simulationResults/relnonlin_div+envvar_cv_unstable2stable.RDS"
rnonlin_unstable2stable_file <- "../simulationResults/relnonlin_div+envvar_cv_unstable2stable_collated.RDS"
if(file.exists(rnonlin_unstable2stable_file) == FALSE) {
  collate_sims(in_file = rnonlin_unstable2stable_large_file, 
               parameter_matrix = parameter_matrix, 
               save_file = rnonlin_unstable2stable_file)
}



##  Read in collated results and do some housekeeping columns
rnonlin_stable2unstable_cvs <- readRDS(rnonlin_stable2unstable_file) 
rnonlin_stable2unstable_cvs$spporder <- "stable2unstable"
rnonlin_unstable2stable_cvs <- readRDS(rnonlin_unstable2stable_file)
rnonlin_unstable2stable_cvs$spporder <- "unstable2stable"
rnonlin_cvs <- rbind(rnonlin_stable2unstable_cvs, rnonlin_unstable2stable_cvs)

##  Calculate initial species richness for each simulation
rnonlin_cvs$init_spprich <- rowSums(rnonlin_cvs[c("N1","N2","N3","N4")])
ids_to_keep <- which(rnonlin_cvs$spprichness == rnonlin_cvs$init_spprich)
rnonlin_cvs_eqrich <- rnonlin_cvs[ids_to_keep,]

##  Summarise by taking mean CV at different realized richness
rnonlin_cvs_mean_full <- ddply(rnonlin_cvs_eqrich, .(Rsd_annual, spprichness, spporder), 
                               summarise,
                               cv = mean(cv))
rnonlin_cvs_mean <- subset(rnonlin_cvs_mean_full, Rsd_annual<1.2 & Rsd_annual>0  & spprichness>0)

##  Plot and save
## Find species first occurences for lines
firstones_rnonlin <- ddply(rnonlin_cvs_mean, 
                           .(spporder, spprichness), 
                           summarise, 
                           Rsd=min(Rsd_annual),
                           cv=cv[which.min(Rsd_annual)])

##  Plot and save
ggplot()+
  geom_segment(data=firstones_rnonlin, aes(x=Rsd, y=cv, xend=Rsd, yend=0, color=as.factor(spprichness)), size=0.5, alpha=0.7)+
  geom_line(data=rnonlin_cvs_mean, aes(x=Rsd_annual, y=cv, color=as.factor(spprichness)), size=0.7, alpha=0.7)+
  geom_point(data=rnonlin_cvs_mean, aes(x=Rsd_annual, y=cv, color=as.factor(spprichness)), size=2, alpha=1)+
  geom_point(data=rnonlin_cvs_mean, aes(x=Rsd_annual, y=cv, color=as.factor(spprichness)), size=2, alpha=1, shape=1, color="grey15")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  facet_wrap("spporder", ncol=2)+
  theme_bw()+
  my_theme+
  theme(legend.position = c(0.1, 0.7))
ggsave(filename = "../manuscript/components/relative_nonlinearity_div+envar.png", width = 110, height=60, units = "mm", dpi = 600)

## Log-log plot and save
ggplot(rnonlin_cvs_mean, aes(x=log(Rsd_annual), y=log(cv), color=as.factor(spprichness)))+
  stat_smooth(method="lm", size=0.7, alpha=0.7, se=FALSE)+
  geom_point(size=2, alpha=1)+
  geom_point(size=2, alpha=1, shape=1, color="grey15")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("SD of annual resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  facet_wrap("spporder", ncol=2)+
  theme_bw()+
  my_theme+
  theme(legend.position = c(0.95, 0.15))
ggsave(filename = "../manuscript/components/relative_nonlinearity_div+envar_loglog.png", width = 8, height=4, units = "in", dpi = 82)

##  Get slopes and plot
rnonlin_cvs_slopes <- list()
for(do_order in unique(rnonlin_cvs_mean$spporder)){
  for(do_spprich in unique(rnonlin_cvs_mean$spprichness)){
    tmp_df <- subset(rnonlin_cvs_mean, spporder==do_order & spprichness==do_spprich)
    if(nrow(tmp_df)>2){
      tmp_slope <- get_slope_rnonlin(tmp_df)
      tmpout <- data.frame(spporder=do_order, spprich=do_spprich, 
                           slope=tmp_slope[1], lowslope=tmp_slope[2], hislope=tmp_slope[3])
      rnonlin_cvs_slopes <- rbind(rnonlin_cvs_slopes, tmpout)
    }
  }
}

ylabel <- bquote("Slope (CV vs." ~ sigma[R] ~ ")")
ggplot(rnonlin_cvs_slopes, aes(x=spprich, y=slope, color=as.factor(spporder)))+
  geom_ribbon(aes(ymin=lowslope, ymax=hislope, fill=as.factor(spporder)), alpha=0.2, color=NA)+
  geom_point()+
  geom_line()+
  ylab(ylabel)+
  xlab("Realized species richness")+
  scale_fill_manual(values=c("red", "blue"), name="Spp. Addition Order", labels=c("Stable to unstable", "Unstable to stable"))+
  scale_color_manual(values=c("red", "blue"), name="Spp. Addition Order", labels=c("Stable to unstable", "Unstable to stable"))+
  theme_bw()+
  my_theme
ggsave("../manuscript/components/relative_nonlinearity_div+envar_loglog_slopes.png", width = 4, height=2, units = "in", dpi = 120)

