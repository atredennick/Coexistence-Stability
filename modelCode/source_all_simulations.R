##	Source files for Tredennick et al.
##	"The mechanism of coexistence shapes the diversity-stability relationship"
##	
##	These scripts reproduce all our results, and are identified
##	by figure number below. Each script will take an hour or more
##	to run. So proceed with caution.
##
##	Many of these scripts generate model output that will be saved as .RDS
##	files, but they do not plot the results. See our code archive
##	(link in the paper) for our plotting code. Plotting code includes
##	code for generating SI material plots.



####
####	FIGURE 2 -- Diversity-Stability Relationships
####
##	Figure 2A
source("storageeffect_sims_multispecies.R")

##	Figure 2B
source("storageeffect_sims_div_stability.R")

##	Figure 2C
source("relativenonlin_sims_multispecies.R")

##	Figure 2D
source("relnonlin_sims_div_stability.R")



####
####	FIGURE 3 -- Effect of Environmental Variability (Storage Effect)
####
source("storageeffect_sims_div+envvar_stability_varycomp.R")



####
####	FIGURE 4 -- Effect of Environmental Variability (Relative Nonlinearity)
####
##	Figure 4A
source("relnonlin_sims_div+envvar_stability.R")

##	Figure 4B
source("relnonlin_sims_div+envvar_stability_revpool.R")



####
####	FIGURE 5 -- Example Time Series
####
source("infographic_sims.R")



