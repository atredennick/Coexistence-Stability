################################################################################
##  relative_nonlinearity_fourspecies_Rstars.R: makes a plot of four species' ##
##  resource uptake curves when coexisting via relative nonlinearity          ##
################################################################################


# Growth function parameters
grow_parameters <- list (
  r = c(1,5,10,25),           # max growth rate for each species
  a = c(2,5,10,25),           # rate parameter for Hill function 
  b = c(2.5,20,30,45),   # shape parameter for Hill function
  eps = c(0.5,0.5,0.5,0.5)  # resource-to-biomass efficiency
)

##  Resource uptake function (Hill function)
uptake_R <- function(r, R, a, b) {
  return((r*R^a) / (b^a + R^a))
}

get_Ruptakes <- function(parms){
  R <- seq(0,100,1)
  out_r <- matrix(ncol=4, nrow=length(R))
  alpha <- parms$a
  beta <- parms$b
  for(i in 1:nrow(out_r)){
    out_r[i,1] <- uptake_R(parms$r[1], R[i], alpha[1], beta[1])
    out_r[i,2] <- uptake_R(parms$r[2], R[i], alpha[2], beta[2])
    out_r[i,3] <- uptake_R(parms$r[3], R[i], alpha[3], beta[3])
    out_r[i,4] <- uptake_R(parms$r[4], R[i], alpha[4], beta[4])
  }
  uptake <- data.frame(species=rep(c(1:4), each=nrow(out_r)),
                       resource=rep(R, times=4),
                       uptake=c(out_r[,1], out_r[,2], out_r[,3], out_r[,4]))
  return(uptake)
}

outs <- get_Ruptakes(grow_parameters)
mycols <- c("#009953",
            "#a55bcf",
            "#efb332",
            "#4a005f")
ggplot(outs, aes(x=resource, y=uptake, color=as.factor(species)))+
  geom_line(size=1)+
  scale_color_manual(values=mycols)+
  theme_few()+
  guides(color=FALSE)+
  xlab("Resource Level")+
  ylab("Resource Uptake")
ggsave("../manuscript/components/fourspp_Ruptake_relnonlin.png", width = 4, height = 4, units = "in", dpi = 100)

