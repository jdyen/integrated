# matrix model modules

# stage based
stage <- function(nstage, nsite) {
  
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
age <- function(nstage, nsite) {
  
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
unstructured <- function(nstage, nsite) {
  
  # hyperpriors for sds
  demo_sd <- lognormal(mean = 0.0, sd = 3.0, dim = (nstage + 1))
  
  # fecundity and survival priors
  params <- vector("list", length = (nstage * nstage))
  params[[1]] <- ilogit(normal(mean = 0.0, sd = demo_sd[1], dim = nstage))
  for (i in 2:nstage) {
    params[[i]] <- lognormal(mean = 0.0, sd = demo_sd[2], dim = nstage)
  }
  for (i in (nstage + 1):length(params)) {
    params[[i]] <- ilogit(normal(mean = 0.0,
                                 sd = demo_sd[(floor((i - 1) / nstage) + 2)],
                                 dim = nstage))
  }
  
  # collate outputs
  out <- list(params = params,
              demo_sd = demo_sd)
  
  # return outputs
  out
  
}
