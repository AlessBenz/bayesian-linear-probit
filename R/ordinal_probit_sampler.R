# Ordinal probit Gibbs sampler with MH updates for cutpoints
# Dependencies: truncnorm, MASS

update_z <- function(y, X, beta, gamma) {
  n <- length(y)
  eta <- drop(X %*% beta)
  z <- numeric(n)
  
  for (i in seq_len(n)) {
    z[i] <- truncnorm::rtruncnorm(
      1, mean = eta[i], sd = 1,
      a = gamma[y[i] + 1],
      b = gamma[y[i] + 2]
    )
  }
  z
}

update_beta <- function(z, X, beta_prior_mean, beta_prior_cov) {
  XtX <- crossprod(X)
  XtZ <- crossprod(X, z)
  
  V0_inv <- solve(beta_prior_cov)
  V_post <- solve(XtX + V0_inv)
  m_post <- V_post %*% (XtZ + V0_inv %*% beta_prior_mean)
  
  drop(MASS::mvrnorm(1, mu = drop(m_post), Sigma = V_post))
}

loglik_y_given_beta_gamma <- function(gamma, y, X, beta) {
  eta <- drop(X %*% beta)
  upper <- gamma[y + 2] - eta
  lower <- gamma[y + 1] - eta
  diff <- stats::pnorm(upper) - stats::pnorm(lower)
  sum(log(pmax(diff, .Machine$double.eps)))
}

proposal_gamma <- function(gamma, sigma) {
  g <- gamma
  L <- length(g)
  for (j in 3:(L - 1)) {  # updates gamma_2 ... gamma_5
    g[j] <- truncnorm::rtruncnorm(
      1, mean = gamma[j], sd = sigma,
      a = g[j - 1],
      b = gamma[j + 1]
    )
  }
  g
}

log_q_gamma <- function(g_to, g_from, sigma) {
  L <- length(g_from)
  logq <- 0
  for (j in 3:(L - 1)) {
    a <- g_to[j - 1]
    b <- g_from[j + 1]
    dens <- truncnorm::dtruncnorm(g_to[j], mean = g_from[j], sd = sigma, a = a, b = b)
    logq <- logq + log(pmax(dens, .Machine$double.eps))
  }
  logq
}


update_gamma <- function(gamma, beta, y, X, sigma_MH) {
  g_new <- proposal_gamma(gamma, sigma_MH)
  
  logR <- loglik_y_given_beta_gamma(g_new, y, X, beta) -
    loglik_y_given_beta_gamma(gamma, y, X, beta) +
    log_q_gamma(gamma, g_new, sigma_MH) -
    log_q_gamma(g_new, gamma, sigma_MH)
  
  if (log(stats::runif(1)) < logR) gamma <- g_new
  
  gamma[1] <- -Inf
  gamma[2] <- 0
  gamma[length(gamma)] <- Inf
  gamma
}

gibbs_mh <- function(y, X, init_beta, init_gamma, n_iter,
                     beta_prior_mean, beta_prior_cov, sigma_gamma) {
  beta <- init_beta
  gamma <- init_gamma
  n <- length(y)
  
  samples_beta <- matrix(NA_real_, n_iter, length(beta))
  samples_gamma <- matrix(NA_real_, n_iter, length(gamma))
  
  z <- numeric(n)
  
  for (iter in seq_len(n_iter)) {
    gamma <- update_gamma(gamma, beta, y, X, sigma_gamma)
    z <- update_z(y, X, beta, gamma)
    beta <- update_beta(z, X, beta_prior_mean, beta_prior_cov)
    
    samples_beta[iter, ] <- beta
    samples_gamma[iter, ] <- gamma
  }
  
  list(beta = samples_beta, gamma = samples_gamma)
}

