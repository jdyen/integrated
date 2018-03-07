# full greta model fit
estimate_mpm <- function(pop_samp,
                         site,
                         growdat = NULL,
                         mat_type = "stage",
                         dens_depend = "none",
                         greta_settings = list()) {
  
  # work out number of sites
  nsite <- length(unique(site))
  
  # calculate the number of stages
  nstage <- nrow(pop_samp[[1]])
  
  # extract parameters for model fit
  nsamp <- sapply(pop_samp, ncol)
  
  # set growth data
  growdat_trans <- otolith_prior_fun(growdat, nstage)
  
  # define greta parameters
  mat_params <- get(mat_type)(nstage = nstage, nsite = nsite, growdat = growdat)
  params <- mat_params$params
  demo_sd <- mat_params$demo_sd
  
  # mean initial pop size list
  mu_est_init <- vector("list", length = length(pop_samp))
  
  # mean pop size list
  mu_est <- vector("list", length = length(pop_samp))
  
  # define mpm matrix
  mpm_mat <- vector("list", length = nsite)
  
  # setup a parameter for the density dependence functions (ignored if dens_depend = "none")
  dens_param <- lognormal(mean = 0.0, sd = 1.0, dim = 1)
  
  # loop over each sampled pop and define full dynamic model
  for (j in seq_len(nsite)) {
    mpm_mat[[j]] <- array(do.call("c", lapply(params, function(x) x[j])),
                          dim = c(nstage, nstage))
    for (k in seq_len(ncol(growdat_trans))) {
      
      # set likelihood on growth data
      growdat_tmp <- matrix(growdat_trans[, k], ncol = nstage)
      distribution(growdat_tmp) <- multinomial(size = sum(growdat_tmp),
                                               prob = mpm_mat[[j]][, k],
                                               dim = 1)
      
    }
  }
  for (j in seq_along(pop_samp)) {
    mu_est_init[[j]] <- lognormal(meanlog = 0.0, sdlog = 2.0, dim = nstage)
    mu_est[[j]] <- iterate_state(t(mpm_mat[[site[j]]]),
                                 mu_est_init[[j]],
                                 dens_param,
                                 seq_len(nsamp[j]),
                                 dens_form = dens_depend)
  }
  
  # likelihood
  mu_init_vec <- do.call("c", mu_est_init)
  param_vec <- do.call("c", params)
  y_vec <- unlist(pop_samp)
  mu_vec <- do.call("c", mu_est)
  distribution(y_vec) <- poisson(mu_vec)
  
  # define model
  mod <- model(mu_vec,
               mu_init_vec,
               param_vec,
               demo_sd,
               dens_param)
  
  # set greta settings
  greta_set <- list(nsamples = 1000,
                    nwarmup = 1000,
                    inits = "random")
  greta_set[names(greta_settings)] <- greta_settings
  greta_set$inits <- switch(greta_set$inits,
                            "random" = rnorm(length(mod$dag$example_parameters())),
                            rep(0, length(mod$dag$example_parameters())))
  
  # sample from model
  samples <- mcmc(mod,
                  n_samples = greta_set$nsamples,
                  warmup = greta_set$nwarmup,
                  initial_values = greta_set$inits)
  
  # return outputs
  out <- list(samples = samples,
              model = mod,
              pop_samp = pop_samp,
              site = site,
              mat_type = mat_type,
              dens_depend = dens_depend)
  out
  
}
