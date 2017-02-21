##	Source files for Tredennick, Adler, & Adler
##	  "The relationship between species richness and ecosystem variability
##     is shaped by the mechanism of coexistence."
##	
##	These scripts reproduce all our results, and are identified
##	by figure number below. Each script will take an hour or more
##	to run. So proceed with caution. You may also run into file path
##  problems if you don't have the appropriate directory structure. But,
##  if you unzipped this entire fileset, you should be good to go by
##  setting your working directory to this source file location.
##
##	Many of these scripts generate model output that will be saved as .RDS
##	files, but they do not plot the results. Plotting code includes
##	code for generating SI material plots. Most of the figures required
##  post-production work in Adobe or similar to generate the figures as they
##  appear in the published version of the paper.
##
##  GITHUB REPO: http://github.com/atredennick/Coexistence-Stability
##
##  CONTACT:     Andrew Tredennick (atredenn@gmail.com)
##  LAST UPDATE: Feb. 21, 2017



####
####	FIGURE 2 -- Diversity-Stability Relationships ----
####
##	Figure 2A simulation script
source("storageeffect_sims_multispecies.R")

##	Figure 2B simulation script
source("storageeffect_sims_div_stability.R")

##	Figure simulation script
source("relativenonlin_sims_multispecies.R")

##	Figure simulation script
source("relnonlin_sims_div_stability.R")

##  Generate Figure 2 Panels
source("plot_richness_variability.R")



####
####	FIGURE 3 -- Effect of Environmental Variability (Storage Effect) ----
####
source("storageeffect_sims_div+envvar_stability_varycomp.R")



####
####	FIGURE 4 -- Effect of Environmental Variability (Relative Nonlinearity) ----
####
##	Figure 4A
source("relnonlin_sims_div+envvar_stability.R")

##	Figure 4B
source("relnonlin_sims_div+envvar_stability_revpool.R")



####
####	FIGURE 5 -- Example Time Series ----
####
source("infographic_sims.R")



