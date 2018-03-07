# construct fitted matrices from model outputs
construct_fitted <- function(params, nstage, nsite) {
  
  out <- array(NA, dim = c(nstage, nstage, nsite))
  
  out_tmp <- array(params, dim = c(nsite, nstage, nstage))
  for (j in seq_len(nsite)) {
    out[, , j] <- t(out_tmp[j, , ])
  }
  
  # return outputs
  out
  
}
