# integrated model with data on individual growth and population
#   size distributions

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

# fit model with correct density dependence form
mod <- estimate_mpm(pop_samp = pop_data,
                    site = rep(1, length(pop_data)),
                    growdat = growth_data,
                    mat_type = "stage",
                    dens_depend = "none",
                    greta_settings = list(nsamples = 100,
                                          nwarmup = 100,
                                          inits = "random"))

# test other mat_types and dens_depend
