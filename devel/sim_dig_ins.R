##  Development script to dig into simulations

rm(list=ls())

path2sims <- "../simulationResults/"
seasons_to_exclude <- 500

# Look at the storage effect and relative nonlinearity combo results
both_varvars <- matrix(c(-0.5,1,2,2.5,5,5,20,
                         -0.5,5,5,20,5,5,20,
                         0.5,1,2,2.5,5,5,20,
                         0.5,5,5,20,5,5,20), ncol=7, byrow = TRUE)
colnames(both_varvars) <- c("rho", "r1", "a1", "b1", "r2", "a2", "b2")

both_sims <- readRDS(paste0(path2sims,"storageeffect_relativenonlinearity_equilibrium_runs.RDS"))
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
both_community_biomass_cv <- ddply(save_seasons_both, .(rel_nonlin, rho), summarise,
                                   cv_biomass = sd(total_biomass)/mean(total_biomass),
                                   cv_n1 = sd(N1)/mean(N1),
                                   cv_n2 = sd(N2)/mean(N2),
                                   mu_n1 = mean(N1),
                                   mu_n2 = mean(N2),
                                   sd_n1 = sd(N1),
                                   sd_n2 = sd(N2),
                                   mu_total = mean(total_biomass),
                                   sd_total = sd(total_biomass))

### Total biomass higher when species negatively correlated due to competitive release of N1 during high resource years
### and standard deviation lower b/c comp not causing fluctuations??? Loreau paper???

full_melt <- melt(save_seasons_both, id.vars = c("rho", "rel_nonlin", "simnum", "timestep"))
ggplot(subset(full_melt, variable %in% c("N1", "N2") & timestep%in%c(500:550)))+
  geom_line(aes(x=timestep, y=value, color=variable))+
  facet_grid(rho~rel_nonlin)
