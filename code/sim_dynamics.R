# updaters for pop dynamics

# no density dependence
none <- function(nsim, pop, popmat, dens_param = NULL) {
  
  for (k in 2:nsim) {
    pop[, k] <- rpois(n = nrow(pop),
                      lambda = crossprod(t(popmat),
                                         pop[, (k - 1)]))
  }
  
  # return outputs
  pop
  
}

# Beverton-Holt
bh <- function(nsim, pop, popmat, dens_param = NULL) {
  
  dens_param <- ifelse(is.null(dens_param), 1, dens_param)
  
  for (k in 2:nsim) {
    nkm1 <- sum(as.numeric(pop[, (k - 1)]))
    scale_factor <- nkm1 / (1 + dens_param * nkm1)
    pop[, k] <- rpois(n = nrow(pop),
                      lambda = crossprod((scale_factor * t(popmat)),
                                         pop[, (k - 1)]))
  }
  
  # return outputs
  pop
  
}

# Ricker
ricker <- function(nsim, pop, popmat, dens_param = NULL) {
  
  dens_param <- ifelse(is.null(dens_param), 1, dens_param)
  
  for (k in 2:nsim) {
    nkm1 <- sum(pop[, (k - 1)])
    scale_factor <- nkm1 * exp(-dens_param * nkm1)
    scale_factor <- ifelse(scale_factor < 0.1, 0.1, scale_factor)
    pop[, k] <- rpois(n = nrow(pop),
                      lambda = crossprod((scale_factor * t(popmat)),
                                         pop[, (k - 1)]))
  }
  
  # return outputs
  pop
  
}
