summarise_mpm <- function(fitted) { 
  
  # unpack fitted list
  samples <- fitted$samples
  pop_samp <- fitted$pop_samp
  sites <- fitted$site
  mat_type <- fitted$mat_type
  dens_depend <- fitted$dens_depend
  nsamp <- sapply(pop_samp, ncol)
  
  # work out length of pollution gradient
  nsite <- length(unique(sites))
  
  # calculate the number of stages
  nstage <- nrow(pop_samp[[1]])
  
  # calculate mean parameter estimates from the posterior distribution
  params <- apply(samples[[1]], 2, mean)
  params_lower <- apply(samples[[1]], 2, quantile, 0.1)
  params_upper <- apply(samples[[1]], 2, quantile, 0.9)
  
  # extract the survival and fecundity parameters
  demo_params <- params[grep("param_vec", names(params))]
  demo_params_lower <- params_lower[grep("param_vec", names(params))]
  demo_params_upper <- params_upper[grep("param_vec", names(params))]

  mat_est <- construct_fitted(params = demo_params,
                              nstage = nstage,
                              nsite = nsite)
  mat_est_lower <- construct_fitted(params = demo_params_lower,
                                    nstage = nstage,
                                    nsite = nsite)
  mat_est_upper <- construct_fitted(params = demo_params_upper,
                                    nstage = nstage,
                                    nsite = nsite) 
  
  # extract fitted values
  fitted_mean <- params[grep("mu_vec", names(params))]
  fitted_lower <- params_lower[grep("mu_vec", names(params_lower))]
  fitted_upper <- params_upper[grep("mu_vec", names(params_upper))]
  
  # reformat fitted values
  fitted_mean <- lapply(seq_along(pop_samp),
                        function(i) matrix(fitted_mean[((i - 1) * (nstage * nsamp[i]) + 1):(i * nstage * nsamp[i])],
                                           ncol = nsamp[i]))
  fitted_lower <- lapply(seq_along(pop_samp),
                         function(i) matrix(fitted_lower[((i - 1) * (nstage * nsamp[i]) + 1):(i * nstage * nsamp[i])],
                                            ncol = nsamp[i]))
  fitted_upper <- lapply(seq_along(pop_samp),
                         function(i) matrix(fitted_upper[((i - 1) * (nstage * nsamp[i]) + 1):(i * nstage * nsamp[i])],
                                            ncol = nsamp[i]))
  
  # return everything
  out <- list(mat_est = mat_est,
              mat_est_lower = mat_est_lower,
              mat_est_upper = mat_est_upper,
              fitted_mean = fitted_mean,
              fitted_lower = fitted_lower,
              fitted_upper = fitted_upper)
  
  # add density dependence parameter
  if (fitted$dens_depend != "none")
    out$dens_param <- params[grep("dens_param", names(params))]
  
  out
  
}

project_fitted_mpm <- function(mod, niter = 10) {
  
  mod_summary <- summarise_mpm(mod)
  nsite <- length(mod_summary$fitted_mean)
  nstage <- nrow(mod_summary$mat_est[, , 1])
  
  dd_param <- NULL
  if (mod$dens_depend != "none")
    dd_param <- mod_summary$dens_param
  
  pop <- sim_pop <- vector("list", length = nsite)
  for (i in seq_len(nsite)) {
    
    abund_fitted <- mod_summary$fitted_mean[[i]]
    
    sim_pop[[i]] <- matrix(NA, nrow = nstage, ncol = (niter + ncol(abund_fitted)))
    pop[[i]] <- matrix(NA, nrow = nstage, ncol = niter)
    pop[[i]][, 1] <- rpois(n = nstage, lambda = abund_fitted[, ncol(abund_fitted)])
    pop[[i]] <- get(mod$dens_depend)(nsim = niter, pop = pop[[i]],
                                     popmat = mod_summary$mat_est[, , mod$site[i]],
                                     dens_param = dd_param)  
    sim_pop[[i]] <- cbind(abund_fitted, pop[[i]])
    
  } 
  
  sim_pop
  
}

plot_fitted_mpm <- function(mod)  {
  
  # extract necessary plot values
  mod_summary <- summarise_mpm(mod)
  pop_samp <- mod$pop_samp
  
  # key stats
  nstage <- nrow(pop_samp[[1]])
  nsamp <- sapply(pop_samp, ncol)
  site <- mod$site
  
  # fitted vals
  fitted_mean <- mod_summary$fitted_mean
  fitted_upper <- mod_summary$fitted_upper
  fitted_lower <- mod_summary$fitted_lower
  
  # set up plot parameters
  old_par <- list(mar = par()$mar, mfrow = par()$mfrow)
  par(mfrow = c(3, 1), mar = c(5, 5, 2.5, 1))
  
  # subset for stages
  plot_mean <- vector("list", length = nstage)
  plot_upper <- vector("list", length = nstage)
  plot_lower <- vector("list", length = nstage)
  plot_real <- vector("list", length = nstage)
  
  # set colours
  col_pal_main <- viridis::inferno(256, alpha = 1)[c(20, seq(70, 200,
                                                             length = length(mod$site)))]
  col_pal_sub <- viridis::inferno(256, alpha = 0.4)[c(20, seq(70, 200,
                                                              length = length(mod$site)))]
  
  # plot
  plot_lab <- paste("Stage", 1:nstage, sep = " ")
  if (nstage == 3)
    plot_lab <- c("Larvae", "Juveniles", "Adults")
  for (i in seq_len(nstage)) {
    
    plot_mean <- lapply(fitted_mean, function(x) x[i, ])
    plot_lower <- lapply(fitted_lower, function(x) x[i, ])
    plot_upper <- lapply(fitted_upper, function(x) x[i, ])
    plot_real <- lapply(pop_samp, function(x) x[i, ])
    
    ylims <- range(c(unlist(plot_mean),
                     unlist(plot_real),
                     unlist(plot_lower),
                     unlist(plot_upper)))
    
    plot(pop_samp[[1]][1, ] ~ c(1:max(nsamp)),
         type = "n", bty = "l", las = 1,
         xaxt = "n", yaxt = "n",
         ylim = ylims, xlim = c(1, max(nsamp)),
         xlab = "", ylab = "")
    
    for (j in seq_along(plot_mean)) {
      
      polygon(c(1:nsamp[j], nsamp[j]:1),
              c(plot_upper[[j]], rev(plot_lower[[j]])),
              border = NA, col = col_pal_sub[site[j] + 1])
      lines(plot_real[[j]] ~ c(1:nsamp[j]),
            lwd = 1.2, col = col_pal_main[site[j] + 1])
      
    }
    
    axis(1, at = c(1:max(nsamp)))
    mtext("Sample", side = 1, adj = 0.5, line = 2.5)
    axis(2)
    mtext("Abundance", side = 2, adj = 0.5, line = 3)
    mtext(plot_lab[i], side = 3, adj = 0.01, line = 0.5)
    
  }
  
  par(mar = old_par$mar, mfrow = old_par$mfrow)
  
}

plot_summed_abund <- function(mod) {
  
  # extract necessary plot values
  mod_summary <- summarise_mpm(mod)
  fitted_mean <- mod_summary$fitted_mean[[1]]
  fitted_upper <- mod_summary$fitted_upper[[1]]
  fitted_lower <- mod_summary$fitted_lower[[1]]
  pop_samp <- mod$pop_samp[[1]]
  
  # set colours
  col_pal_main <- viridis::inferno(256, alpha = 1)[c(20, seq(70, 200,
                                                             length = nrow(pop_samp)))]
  col_pal_sub <- viridis::inferno(256, alpha = 0.4)[c(20, seq(70, 200,
                                                              length = nrow(pop_samp)))]
  
  x_vals <- c(1999:2017)
  
  plot_mean <- apply(fitted_mean, 2, sum)
  plot_lower <- apply(fitted_lower, 2, sum)
  plot_upper <- apply(fitted_upper, 2, sum)
  plot_real <- apply(pop_samp, 2, sum)
  
  plot(pop_samp[1, ] ~ x_vals,
       type = "n", bty = "l", las = 1,
       xaxt = "n", yaxt = "n",
       ylim = range(c(plot_mean, plot_lower, plot_upper, plot_real)),
       xlim = range(x_vals),
       xlab = "", ylab = "")
  
  polygon(c(x_vals, rev(x_vals)),
          c(plot_upper, rev(plot_lower)),
          border = NA, col = col_pal_sub[1])
  lines(plot_real ~ x_vals,
        lwd = 1.2, col = col_pal_main[1])
  
  axis(1, at = x_vals[seq(1, length(x_vals), by = 2)], las = 1)
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  axis(2, las = 1)
  mtext("Abundance", side = 2, adj = 0.5, line = 3.5)
  
}

plot_matrix <- function(mod) {
  
  mod_summary <- summarise_mpm(mod)
  image(t(mod_summary$mat_est[, , 1]),
        x = seq(0, 54000, length = ncol(mod_summary$mat_est[, , 1])),
        y = seq(0, 54000, length = ncol(mod_summary$mat_est[, , 1])),
        col = viridis::inferno(256)[c(1, 51:256)],
        breaks = c(seq(0, 1, length = 151), seq(1, 7, length = 57)),
        xlab = "",
        ylab = "",
        las = 1)
  mtext("Size (g) at time t - 1", side = 1, line = 2.5, cex = 1)
  mtext("Size (g) at time t", side = 2, line = 3.8, cex = 1)
  
}

plot_projections <- function(mod, niter = 25, nsim = 30) {

  out <- project_fitted_mpm2(mod, niter = niter)[[1]]
  x_vals <- c(1999:(1998 + ncol(out)))

  mod_summary <- summarise_mpm(mod)
  fitted_mean <- mod_summary$fitted_mean[[1]]
  fitted_upper <- mod_summary$fitted_upper[[1]]
  fitted_lower <- mod_summary$fitted_lower[[1]]
  pop_samp <- mod$pop_samp[[1]]
  
  plot_mean <- apply(out, 2, sum)
  plot_real <- apply(pop_samp, 2, sum)
  plot_lower <- apply(fitted_lower, 2, sum)
  plot_upper <- apply(fitted_upper, 2, sum)
  
  plot(plot_mean ~ x_vals,
       type = "n", bty = "l", las = 1,
       xaxt = "n", yaxt = "n",
       ylim = c(0, 5000),
       xlim = range(x_vals),
       xlab = "", ylab = "")
  
  # set colours
  col_pal_main <- viridis::inferno(256, alpha = 1)[c(20, seq(70, 200,
                                                             length = nsim))]
  col_pal_sub <- viridis::inferno(256, alpha = 0.4)[c(20, seq(70, 200,
                                                              length = nsim))]
  
  for (i in seq_len(nsim)) {
    
    out <- project_fitted_mpm2(mod, niter = niter, dd_param = mod_summary$dens_param)[[1]]
    plot_mean <- apply(out, 2, sum)
    lines(plot_mean[(length(plot_real)):length(plot_mean)] ~ x_vals[(length(plot_real)):length(plot_mean)],
          lwd = 1.2, col = col_pal_main[i])
    
  }
  
  polygon(c(x_vals[seq_along(plot_upper)], rev(x_vals[seq_along(plot_upper)])),
          c(plot_upper, rev(plot_lower)),
          border = NA, col = col_pal_sub[1])
  lines(plot_mean ~ x_vals,
        lwd = 1.2, col = col_pal_main[1])
  lines(plot_real ~ x_vals[seq_along(plot_real)],
        lwd = 1.2, col = col_pal_main[1])
  
  axis(1, at = x_vals[seq(1, length(x_vals), by = 2)], las = 1)
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  axis(2, las = 1)
  mtext("Abundance", side = 2, adj = 0.5, line = 3.5)

}
