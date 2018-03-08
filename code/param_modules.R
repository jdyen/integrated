# matrix model modules

# stage based
stage <- function(nstage, nsite, growdat) {
  
  # hyperpriors for sds
  demo_sd <- lognormal(mean = 0.0,
                       sd = 1.0,
                       dim = 2)
  
  # fecundity and survival priors
  params <- vector("list", length = (nstage * nstage))
  params[[1]] <- ilogit(normal(mean = 0.0,
                               sd = demo_sd[1],
                               dim = nsite))
  for (i in 2:(nstage - 1)) {
    params[[i]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
  }
  params[[nstage]] <- lognormal(mean = 0.0,
                                sd = demo_sd[2],
                                dim = nsite) 
  for (i in 2:nstage) {
    for (j in seq_len(nstage)) {
      if ((j == (i - 1)) | (j == (i - 2)) | (j == i)) {
        params[[((i - 1) * nstage) + j]] <- ilogit(normal(mean = 0.0,
                                                          sd = demo_sd[1],
                                                          dim = nsite))
      } else {
        params[[((i - 1) * nstage) + j]] <- uniform(min = 0.0,
                                                    max = 0.0001,
                                                    dim = nsite)
      }
    }
  }
  
  # collate outputs
  out <- list(params = params,
              demo_sd = demo_sd)
  
  # return outputs
  out
  
}

# age based
age <- function(nstage, nsite, growdat) {
  
  if (!is.null(growdat))
    warning("growdat must be suppllied", call. = FALSE)
  
  # hyperpriors for sds
  demo_sd <- lognormal(mean = 0.0, sd = 3.0, dim = (nstage + 1))
  
  # fecundity and survival priors
  params <- vector("list", length = (nstage * nstage))
  for (i in 1:(nstage - 1)) {
    params[[i]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
  }
  params[[nstage]] <- lognormal(mean = 0.0, sd = demo_sd[1], dim = nsite)
  for (i in 2:nstage) {
    for (j in seq_len(nstage)) {
      if (j == (i - 1)) {
        params[[((i - 1) * nstage) + j]] <- ilogit(normal(mean = 0.0, sd = demo_sd[i],
                                                          dim = nsite))
      } else {
        if ((i == nstage) & (j == nstage)) {
          params[[((i - 1) * nstage) + j]] <- ilogit(normal(mean = 0.0, sd = demo_sd[(nstage + 1)],
                                                            dim = nsite))
        } else {
          params[[((i - 1) * nstage) + j]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
        }
      }
    }
  }
  
  # collate outputs
  out <- list(params = params, demo_sd = demo_sd)
  
  # return outputs
  out
  
}

# unstructured
unstructured <- function(nstage, nsite, growdat) {
  
  # hyperpriors for sds
  demo_sd <- lognormal(mean = 0.0, sd = 3.0, dim = (nstage + 1))
  
  # fecundity and survival priors
  params_fec <- vector("list", length = (nstage * nstage))
  params_surv <- vector("list", length = (nstage * nstage))
  params_disp <- vector("list", length = (nstage * nstage))
  params_surv[[1]] <- ilogit(normal(mean = 0.0, sd = demo_sd[1], dim = nsite))
  params_disp[[1]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
  params_fec[[1]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
  for (i in 2:nstage) {
    params_surv[[i]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
    params_fec[[i]] <- lognormal(mean = 0.0,
                                 sd = demo_sd[((i - 1) %% nstage) + 1],
                                 dim = nsite)
    params_disp[[i]] <- lognormal(mean = 0.0,
                                  sd = demo_sd[nstage + 1],
                                  dim = nsite)
  } 
  
  for (i in (nstage + 1):length(params_surv)) {
    params_fec[[i]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
    params_surv[[i]] <- ilogit(normal(mean = 0.0,
                                      sd = demo_sd[((i - 1) %% nstage) + 1],
                                      dim = nsite))
    if (i %in% seq((nstage + 2), length(params_surv), by = (nstage + 1))) {
      params_disp[[i]] <- lognormal(mean = 0.0,
                                    sd = demo_sd[nstage + 1],
                                    dim = nsite)
    } else { 
      params_disp[[i]] <- uniform(min = 0.0, max = 0.0001, dim = nsite)
    }
  }
  
  # collate outputs
  out <- list(params_surv = params_surv,
              params_fec = params_fec,
              params_disp = params_disp,
              demo_sd = demo_sd)
  
  # return outputs
  out
  
}
