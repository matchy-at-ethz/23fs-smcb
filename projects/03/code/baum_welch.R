
forward_prob <- function() {
  # placeholder
}

backward_prob <0- function() {
  # placeholder
}

baum_welch <- function(Tr, E,
                       state_label, symbol_label, threshold=1e-5, max_iter=1000) {
  
  # get the number of states and symbols
  K <- length(state_label)
  M <- length(symbol_label)
  
  # set Tr and E to all zeros
  Tr <- matrix(0, K, K)
  E  <- matrix(0, K, M)
  
  ll_hid <- -Inf
  ll_obs <- -Inf
  diff <- Inf
  iter <- 0
  while (diff >= threshold && iter < max_iter) {
    
    ll_obs_old <- ll_obs
    # E-step
    # compute the forward and backward probabilities
    
    
    # M-step
    # compute the new transition and emission matrices
    
    # update the log-likelihood
    log_likelihood <- log_likelihood_new
    
    # update the iteration
    iter <- iter + 1
  }
}
