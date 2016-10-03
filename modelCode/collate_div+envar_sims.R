##  collate_div+envvar_sims.R

rm(list=ls(all.names = TRUE))


####
####  Libraries
####
library("plyr")
library("reshape2")
library("ggplot2")
library("ggthemes")
library("viridis")


####
####  Function for calculating CV from list elements
####
get_cv <- function(x){
  seasons_to_exclude <- 500
  x <- as.data.frame(x)
  names(x) <- c("D1", "D2", "D3", "D4", "N1", "N2", "N3", "N4", "R")
  livestates <- grep("N", colnames(x))
  tmp_totbiomass <- rowSums(x[seasons_to_exclude:nrow(x),livestates])
  tmp_cv <- sd(tmp_totbiomass) / mean(tmp_totbiomass)
  tmp_sppavg <- colMeans(x[seasons_to_exclude:nrow(x),livestates])
  tmp_spprich <- length(which(tmp_sppavg > 1))
  
  tmp_out <- c(tmp_cv, tmp_spprich)
}


#######                            #######
####### STORAGE EFFECT SIMULATIONS #######
#######                            #######


####
####  Re-define the parameter matrix 
####  (from storageeffect_sims_div+envar_stability.R)
####
## Define vectors of parameters to vary -- here, initial conditions to vary richness
## Define vectors of parameters to vary -- here, initial conditions to vary richness
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



####
####  Read in the simulation results -- a BIG list
####
# file_to_get <- "../../../Desktop/storageeffect_div+envvar_stability_asymcomp1.RDS"
# # file_to_get <- "../../../Desktop/storageeffect_div+envvar_stability.RDS"
# sims_list <- readRDS(file_to_get)
# list_elements <- length(sims_list)
# 
# 
# 
# ####
# ####  Apply the get_cv function over the simulation list
# ####
# out_cvs <- lapply(sims_list, get_cv) # apply the function
# out_cvs_df <- do.call(rbind.data.frame, out_cvs) # convert to data frame
# colnames(out_cvs_df) <- c("cv", "spprichness") # give the columns names
# cvs_params <- cbind(parameter_matrix, out_cvs_df) # combine with parameter values
# # saveRDS(cvs_params, file = "../simulationResults/storageeffect_div+envvar_cv_collated.RDS")
# saveRDS(cvs_params, file = "../simulationResults/storageeffect_div+envvar_cv_collated_asymmetric.RDS")

####  READ IN COLLATED RESULTS
cvs_params_symm <- readRDS(file = "../simulationResults/storageeffect_div+envvar_cv_collated.RDS")
cvs_params_asymm <- readRDS(file = "../simulationResults/storageeffect_div+envvar_cv_collated_asymmetric.RDS")
cvs_params_symm$comp <- "1symmetrical"
cvs_params_asymm$comp <- "2asymmetrical"
cvs_params <- rbind(cvs_params_symm, cvs_params_asymm)

####
#### Summarize and plot results
####
get_colors <- function(n_cols){
  return(scales::seq_gradient_pal("coral", "dodgerblue", "Lab")(seq(0,1,length.out=n_cols)))
}
mycols <- get_colors(4)

cvs_plots <- ddply(cvs_params, .(sigE, rho, spprichness, comp), summarise,
              cv = mean(cv))

cvs_plots2 <- subset(cvs_plots, sigE!=0)
ggplot(cvs_plots2, aes(x=sigE, y=cv, color=as.factor(spprichness)))+
  #   geom_line(size=0.7)+
  geom_point(size=2, alpha=0.3)+
  geom_point(size=2, alpha=0.5, shape=1)+
  stat_smooth(method="loess", se=FALSE, size=1)+
  # scale_color_manual(values = mycols, name = "Species\nRichness")+
  scale_color_viridis(discrete = TRUE, name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(comp~rho)+
  theme_few()+
  theme(legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.background = element_rect(colour = "white", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(0.95, 0.15))+
  theme(strip.background = element_blank())
ggsave(filename = "../manuscript/components/storage_effect_div+envar.png", width = 8, height=4, units = "in", dpi = 82)




####  
####  PLOT ASSYMETRIC COMP STORAGE EFFECT; VARY ETAS
####
cvs_plots <- ddply(cvs_params_asymm, .(sigE, rho, spprichness,eta1), summarise,
                   cv = mean(cv))
ggplot(cvs_plots, aes(x=sigE, y=cv, color=as.factor(spprichness)))+
  geom_line(size=0.7)+
  geom_point(size=2, alpha=0.3)+
  geom_point(size=2, alpha=0.5, shape=1)+
  # stat_smooth(method="loess", se=FALSE, size=1)+
  # scale_color_manual(values = mycols, name = "Species\nRichness")+
  # scale_color_viridis(discrete = TRUE, name = "Species\nRichness")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(eta1~rho, scales="free_y")+
  theme_bw()+
  theme(legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.background = element_rect(colour = "white", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(0.95, 0.35),
        panel.grid.major = element_line(colour = "grey25", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"))+
  theme(strip.background = element_blank())
ggsave("../manuscript/components/SI_storageeffect_varyetas.png", width = 8, height=4, units = "in", dpi = 82)

## Now log-log to look at slopes
ggplot(cvs_plots, aes(x=log(sigE), y=log(cv), color=as.factor(spprichness)))+
  # geom_line(size=0.7)+
  geom_point(size=2, alpha=0.3)+
  geom_point(size=2, alpha=0.5, shape=1)+
  stat_smooth(method="lm", se=FALSE, size=1)+
  # scale_color_manual(values = mycols, name = "Species\nRichness")+
  # scale_color_viridis(discrete = TRUE, name = "Species\nRichness")+
  scale_color_brewer(palette = "Set1", name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  facet_grid(eta1~rho, scales="free_y")+
  theme_bw()+
  theme(legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.background = element_rect(colour = "white", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(0.95, 0.35),
        panel.grid.major = element_line(colour = "grey25", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"))+
  theme(strip.background = element_blank())
ggsave("../manuscript/components/SI_storageeffect_varyetas_loglog.png", width = 8, height=4, units = "in", dpi = 82)

## Calculate slopes and plot
get_slope <- function(df_subset) {
  mod <- lm(log(cv) ~ log(sigE), data = df_subset)
  slope <- c(coef(mod)[2], confint(mod)[2,1], confint(mod)[2,2])
  return(slope)
}

slope_df <- list()
for(do_eta in unique(cvs_plots$eta1)){
  for(do_rho in c(-1,0)){
    for(do_spprich in unique(cvs_plots$spprichness)){
      tmp_df <- subset(cvs_plots, eta1==do_eta & spprichness==do_spprich & rho==do_rho)
      tmp_slope <- get_slope(tmp_df)
      tmpout <- data.frame(eta=do_eta, spprich=do_spprich, rho=do_rho, 
                           slope=tmp_slope[1], lowslope=tmp_slope[2], hislope=tmp_slope[3])
      slope_df <- rbind(slope_df, tmpout)
    }
  }
}

ylabel <- bquote("Slope (CV vs." ~ sigma[E] ~ ")")

ggplot(slope_df, aes(x=spprich, y=slope, color=as.factor(eta)))+
  geom_ribbon(aes(ymin=lowslope, ymax=hislope, fill=as.factor(eta)), alpha=0.2, color=NA)+
  geom_point()+
  geom_line()+
  facet_wrap("rho")+
  ylab(ylabel)+
  xlab("Realized species richness")+
  scale_fill_manual(values=c("red", "blue"), name=expression(eta))+
  scale_color_manual(values=c("red", "blue"), name=expression(eta))+
  theme_bw()+
  theme(legend.title=element_text(size=12, face="bold"),
        legend.text=element_text(size=10),
        legend.background = element_rect(colour = "white", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_line(colour = "grey25", linetype = "dotted"),
        panel.grid.minor = element_line(colour = "white", linetype = "dotted"))+
  theme(strip.background = element_blank())
ggsave("../manuscript/components/SI_storageeffect_varyetas_loglog_slopes.png", width = 8, height=4, units = "in", dpi = 82)

#######                                   #######
####### RELATIVE NONLINEARITY SIMULATIONS #######
#######                                   #######

##  Recreate parameter matrix
# DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
#              c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
#              c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
#              c(D=c(1,0,0,0),N=c(1,0,0,0),R=20))

DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(0,1,1,1),N=c(0,1,1,1),R=20),
             c(D=c(0,0,1,1),N=c(0,0,1,1),R=20),
             c(D=c(0,0,0,1),N=c(0,0,0,1),R=20))

## Define vectors of parameters to vary
n_rsd <- 50 # Number of seasonal standard deviation levels
rsd_vec <- pretty(seq(0, 1.4, length.out=n_rsd), n_rsd) # Make a pretty vector
prm <- as.data.frame(rsd_vec)
colnames(prm) <- "Rsd_annual"
# Add in variable parameters to form parameter matrix
DNR_repped <- do.call("rbind", replicate(nrow(prm), DNR,  simplify = FALSE))
prm_repped <- do.call("rbind", replicate(nrow(DNR), prm,  simplify = FALSE))
parameter_matrix <- cbind(DNR_repped, prm_repped)


####
####  Read in the simulation results -- a BIG list
####
# file_to_get <- "../../../Desktop/relnonlin_div+envvar_stability_revpool.RDS"
# sims_list <- readRDS(file_to_get)
# list_elements <- length(sims_list)
# 
# 
# 
# ####
# ####  Apply the get_cv function over the simulation list
# ####
# out_cvs <- lapply(sims_list, get_cv) # apply the function
# out_cvs_df <- do.call(rbind.data.frame, out_cvs) # convert to data frame
# colnames(out_cvs_df) <- c("cv", "spprichness") # give the columns names
# cvs_params <- cbind(parameter_matrix, out_cvs_df) # combine with parameter values
# saveRDS(cvs_params, file = "../simulationResults/relnonlin_div+envvar_cv_revpool_collated.RDS")

####  READ IN COLLATED RESULTS
cvs_params_symm <- readRDS(file = "../simulationResults/relnonlin_div+envvar_cv_collated.RDS")
cvs_params_symm_revpool <- readRDS(file = "../simulationResults/relnonlin_div+envvar_cv_revpool_collated.RDS")
cvs_params_symm$pool <- 1
cvs_params_symm_revpool$pool <- 2
cvs_params <- rbind(cvs_params_symm, cvs_params_symm_revpool)

cvs_plots <- ddply(cvs_params, .(Rsd_annual, spprichness, pool), summarise,
                   cv = mean(cv))

cvs_plots2 <- subset(cvs_plots, Rsd_annual<1.2 & spprichness>0)
ggplot(cvs_plots2, aes(x=Rsd_annual, y=cv, color=as.factor(spprichness)))+
  #   geom_line(size=0.7)+
  geom_point(size=2, alpha=0.3)+
  geom_point(size=2, alpha=0.5, shape=1)+
  stat_smooth(method="loess", se=FALSE, size=1)+
  # scale_color_manual(values = mycols, name = "Species\nRichness")+
  scale_color_viridis(discrete = TRUE, name = "Species\nRichness")+
  xlab(expression(paste("SD of resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  facet_wrap("pool", ncol=2)+
  theme_few()+
  theme(legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.background = element_rect(colour = "white", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(0.08, 0.7))+
  theme(strip.background = element_blank())
ggsave(filename = "../manuscript/components/relative_nonlinearity_div+envar.png", width = 6, height=3, units = "in", dpi = 82)


### Look at plots based on initial pool
cvs_params_symm_revpool$sppstart <- rowSums(cvs_params_symm_revpool[,1:4])
# cvs_plots <- ddply(cvs_params_symm_revpool, .(Rsd_annual, sppstart, pool, spprichness), summarise,
#                    cv = mean(cv))

cvs_plots3 <- subset(cvs_params_symm_revpool, Rsd_annual<1.2)

ggplot(cvs_plots3, aes(x=Rsd_annual, y=cv, color=as.factor(sppstart)))+
  #   geom_line(size=0.7)+
  geom_point(size=2, alpha=0.3)+
  geom_point(size=2, alpha=0.5, shape=1)+
  stat_smooth(method="loess", se=FALSE, size=1)+
  # scale_color_manual(values = mycols, name = "Species\nRichness")+
  scale_color_viridis(discrete = TRUE, name = "Species\nRichness")+
  xlab(expression(paste("SD of resource (",sigma[R], ")")))+
  ylab("CV of community biomass")+
  facet_wrap("pool", ncol=1)+
  theme_few()+
  scale_y_continuous(limits=c(0,2))+
  theme(legend.title=element_text(size=8),
        legend.text=element_text(size=8),
        legend.background = element_rect(colour = "white", size=0.5),
        legend.key = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        legend.position = c(0.2, 0.85))+
  theme(strip.background = element_blank())


#######                             #######
####### EXTRAS EXTRAS EXTRAS EXTRAS #######
#######                             #######

####
####  Look at just one nested species pool
####
# cvs_pool <- rbind(subset(cvs_params, N1==1 & N2==1 & N3==1 & N4==1),
#                   subset(cvs_params, N1==1 & N2==1 & N3==1 & N4==0),
#                   subset(cvs_params, N1==1 & N2==1 & N3==0 & N4==0),
#                   subset(cvs_params, N1==1 & N2==0 & N3==0 & N4==0))
# 
# cvs_pool_avg <- ddply(cvs_pool, .(sigE, spprichness), summarise,
#                       avg_cv = mean(cv))
# 
# ggplot(cvs_pool_avg, aes(x=sigE, y=avg_cv, color=as.factor(spprichness)))+
#   geom_line(size=0.8)+
# #   geom_point(size=2, alpha=0.3)+
# #   geom_point(size=2, alpha=0.5, shape=1)+
#   # stat_smooth(method="lm", se=FALSE, size=1)+
#   scale_color_manual(values = mycols, name = "Species\nRichness")+
#   xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
#   ylab("CV of community biomass")+
#   theme_bw()+
#   theme(legend.title=element_text(size=10),
#         legend.background = element_rect(colour = "grey", size=0.5))



####
####  Look at some time series
####
# ts_list <- readRDS("../../../Desktop/storageeffect_div+envvar_stability.RDS")
# two_spp_ids <- which(parameter_matrix[,1] == 1 & parameter_matrix[,2]==1 & parameter_matrix[,3] == 0 & parameter_matrix[,4]==0)
# pdf("two_spp_timeseries.pdf")
# for(i in 1:length(two_spp_ids)){
#   tsnow <- ts_list[[two_spp_ids[i]]]
#   matplot(tsnow[,c(5:6)], type="l")
# }
# dev.off()

