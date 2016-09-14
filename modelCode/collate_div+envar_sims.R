##  collate_div+envvar_sims.R

rm(list=ls(all.names = TRUE))


####
####  Libraries
####
library("plyr")
library("reshape2")
library("ggplot2")
library("ggthemes")



####
####  Re-define the parameter matrix 
####  (from storageeffect_sims_div+envar_stability.R)
####
## Define vectors of parameters to vary -- here, initial conditions to vary richness
DNR <- rbind(c(D=c(1,1,1,1),N=c(1,1,1,1),R=20),
             c(D=c(1,1,1,0),N=c(1,1,1,0),R=20),
             c(D=c(1,1,0,1),N=c(1,1,0,1),R=20),
             c(D=c(1,0,1,1),N=c(1,0,1,1),R=20),
             c(D=c(0,1,1,1),N=c(0,1,1,1),R=20),
             c(D=c(1,1,0,0),N=c(1,1,0,0),R=20),
             c(D=c(1,0,1,0),N=c(1,0,1,0),R=20),
             c(D=c(1,0,0,1),N=c(1,0,0,1),R=20),
             c(D=c(0,1,1,0),N=c(0,1,1,0),R=20),
             c(D=c(0,1,0,1),N=c(0,1,0,1),R=20),
             c(D=c(0,0,1,1),N=c(0,0,1,1),R=20),
             c(D=c(1,0,0,0),N=c(1,0,0,0),R=20),
             c(D=c(0,1,0,0),N=c(0,1,0,0),R=20),
             c(D=c(0,0,1,0),N=c(0,0,1,0),R=20),
             c(D=c(0,0,0,1),N=c(0,0,1,0),R=20))

n_sig_e <- 50 # Number of cue variance levels
sig_e_vec <- pretty(seq(0, 2, length.out=n_sig_e), n_sig_e) # Make a pretty vector
prm <- as.matrix(sig_e_vec)
colnames(prm) <- "sigE"

# Add in variable parameters to form parameter matrix
DNR_repped <- do.call("rbind", replicate(nrow(prm), DNR,  simplify = FALSE))
prm_repped <- do.call("rbind", replicate(nrow(DNR), prm,  simplify = FALSE))
parameter_matrix <- cbind(DNR_repped, prm_repped)



####
####  Read in the simulation results -- a BIG list
####
file_to_get <- "../simulationResults/storageeffect_div+envvar_stability.RDS"
sims_list <- readRDS(file_to_get)
list_elements <- length(sims_list)



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



####
####  Apply the get_cv function over the simulation list
####
out_cvs <- lapply(sims_list, get_cv) # apply the function
out_cvs_df <- do.call(rbind.data.frame, out_cvs) # convert to data frame
colnames(out_cvs_df) <- c("cv", "spprichness") # give the columns names
cvs_params <- cbind(parameter_matrix, out_cvs_df) # combine with parameter values
saveRDS(cvs_params, file = "../simulationResults/storageeffect_div+envvar_cv_collated.RDS")



####
#### Summarize and plot results
####
get_colors <- function(n_cols){
  return(scales::seq_gradient_pal("red", "blue", "Lab")(seq(0,1,length.out=n_cols)))
}
mycols <- get_colors(4)

ggplot(cvs_params, aes(x=sigE, y=cv, color=as.factor(spprichness)))+
  geom_point(size=2, alpha=0.3)+
  geom_point(size=2, alpha=0.5, shape=1)+
  # stat_smooth(method="lm", se=FALSE, size=1)+
  scale_color_manual(values = mycols, name = "Species\nRichness")+
  xlab(expression(paste("Variance of environmental cue (",sigma[E]^2, ")")))+
  ylab("CV of community biomass")+
  theme_bw()+
  theme(legend.title=element_text(size=10),
        legend.background = element_rect(colour = "grey", size=0.5))








