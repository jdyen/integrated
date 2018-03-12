# fit integrated model with data on individual growth and
#   population size distributions

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
source("./code/install_packages.R")

# source helper functions
source("./code/construct_fitted_matrix.R")
source("./code/helpers.R")
source("./code/mpm_est.R")
source("./code/param_modules.R")
source("./code/sim_dynamics.R")
source("./code/otolith_calc.R")

# load simulated population and growth data
pop_data <- get(load("./data/pop_data.R"))
growth_data <- read.csv("./data/growth_data.csv")

# fit integrated model with stage-structured transitions and
#   Beverton-Holt density dependence.
#   
# Note: 1000 iterations could take >30 min
#
mod_stage <- estimate_mpm(pop_samp = pop_data,
                          site = rep(1, length(pop_data)),
                          growdat = growth_data,
                          mat_type = "stage",  # fit a stage-structured model
                          dens_depend = "bh",  # use Beverton-Holt density dependence
                          greta_settings = list(nsamples = 1000,
                                                nwarmup = 1000,
                                                inits = "random"))

# plot observed and fitted population abundance
plot_summed_abund(mod_stage)

# plot matrix of estimated vital rates
plot_matrix(mod_stage)

# plot projected species's abundances
plot_projections(mod_stage)

# fit a model with age-structured transitions
mod_age <- estimate_mpm(pop_samp = pop_data,
                        site = rep(1, length(pop_data)),
                        growdat = growth_data,
                        mat_type = "age",
                        dens_depend = "none",  # do not include density dependence
                        greta_settings = list(nsamples = 1000,
                                              nwarmup = 1000,
                                              inits = "random"))

# fit a model with unstructured transitions
mod_unstruc <- estimate_mpm(pop_samp = pop_data,
                            site = rep(1, length(pop_data)), 
                            growdat = growth_data,
                            mat_type = "unstructured",
                            dens_depend = "ricker",  # use Ricker density dependence
                            greta_settings = list(nsamples = 1000,
                                                  nwarmup = 1000,
                                                  inits = "random"))
