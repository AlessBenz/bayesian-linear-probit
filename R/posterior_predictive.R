# Posterior predictive simulation utilities
# Note: In an ordinal probit model, the latent variable is z = X beta + epsilon, epsilon ~ N(0,1).
# For PPC, simulate epsilon each time, then discretize using cutpoints.

#' Simulate ordinal outcomes from latent mean and cutpoints
#' @param eta latent mean vector (X %*% beta)
#' @param gamma cutpoints vector with -Inf and +Inf bounds (length K+1)
#' @param include_noise if TRUE, adds N(0,1) noise before discretization (recommended)
#' @return integer vector of categories 0:(K-1)
simulate_ordinal <- function(eta, gamma, include_noise = TRUE) {
  eta <- drop(eta)
  z <- if (include_noise) eta + stats::rnorm(length(eta), 0, 1) else eta
  
  # gamma has length K+1 with -Inf and +Inf; categories are 0:(K-1)
  as.integer(findInterval(z, gamma, rightmost.closed = TRUE) - 1L)
}

#' Generate replicated datasets
#' @param X design matrix
#' @param beta_samples matrix (iter x p)
#' @param gamma_samples matrix (iter x (K+1))
#' @param include_noise include N(0,1) noise on the latent scale
#' @param draw_ids optional indices of posterior draws to use (recommended for speed)
#' @return matrix yrep (length(draw_ids) x n)
simulate_yrep <- function(X, beta_samples, gamma_samples, include_noise = TRUE, draw_ids = NULL) {
  if (is.null(draw_ids)) draw_ids <- seq_len(nrow(beta_samples))
  n <- nrow(X)
  yrep <- matrix(NA_integer_, length(draw_ids), n)
  
  for (s in seq_along(draw_ids)) {
    i <- draw_ids[s]
    eta <- X %*% beta_samples[i, ]
    yrep[s, ] <- simulate_ordinal(eta, gamma_samples[i, ], include_noise = include_noise)
  }
  yrep
}
