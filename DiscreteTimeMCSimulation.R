#############################################
#DISCRETE TIME MONTE CARLO SIM
#############################################

set.seed(42)

m       <- 5
T       <- 300
L       <- 3
z_alpha <- qnorm(0.975)
N_sim   <- 500

#############################################
# TRUE PARAMETERS (fixed across simulations)
#############################################

mu_true <- runif(m, 0.1, 0.3)

alpha_true <- matrix(0, m, m)
alpha_true[2, 1] <- 0.8
alpha_true[3, 2] <- 0.7
alpha_true[4, 3] <- 0.6
alpha_true[5, 4] <- 0.5

beta_true  <- 0.8
sigma_true <- 0.5

#############################################
# HELPER FUNCTIONS
#############################################

alpha_index_map <- function(m) {
  idx <- matrix(NA, m, m)
  counter <- 1
  for (i in 1:m) {
    for (j in 1:m) {
      if (i != j) {
        idx[i, j] <- counter
        counter   <- counter + 1
      }
    }
  }
  idx
}

neg_loglik <- function(par, Y, dtheta, m, L) {
  
  mu      <- exp(par[1:m])
  idx_map <- alpha_index_map(m)
  n_alpha <- m * (m - 1)
  
  eta_alpha <- par[(m + 1):(m + n_alpha)]
  beta      <- exp(par[m + n_alpha + 1])
  sigma     <- exp(par[m + n_alpha + 2])
  
  loglik <- 0
  
  for (i in 1:m) {
    for (b in (L + 1):ncol(Y)) {
      
      lambda <- mu[i]
      
      for (j in 1:m) {
        if (j != i) {
          k        <- idx_map[i, j]
          alpha_ij <- exp(eta_alpha[k])
          
          for (ell in 1:L) {
            lambda <- lambda +
              alpha_ij *
              Y[j, b - ell] *
              exp(-beta * ell) *
              exp(-(dtheta[i, j, b, ell]^2) / (2 * sigma^2))
          }
        }
      }
      
      loglik <- loglik +
        dpois(Y[i, b], lambda = lambda, log = TRUE)
    }
  }
  
  return(-loglik)
}

#############################################
# STORAGE FOR COVERAGE INDICATORS
#############################################

# Parameters to track: beta, sigma, alpha[2,1], alpha[3,2], alpha[4,3], alpha[5,4]
param_names <- c("beta", "sigma", "alpha[2,1]", "alpha[3,2]", "alpha[4,3]", "alpha[5,4]")
true_vals   <- c(beta_true, sigma_true,
                 alpha_true[2,1], alpha_true[3,2],
                 alpha_true[4,3], alpha_true[5,4])

covered     <- matrix(NA, nrow = N_sim, ncol = length(param_names),
                      dimnames = list(NULL, param_names))
converged   <- logical(N_sim)

idx_map <- alpha_index_map(m)
n_alpha <- m * (m - 1)

init_par <- c(
  log(rep(0.2, m)),
  log(rep(0.05, n_alpha)),
  log(0.7),
  log(0.6)
)

#############################################
# MONTE CARLO LOOP
#############################################

cat("Running", N_sim, "simulations...\n")

for (sim in 1:N_sim) {
  
  if (sim %% 50 == 0) cat("  Simulation", sim, "/", N_sim, "\n")
  
  ## ---- Simulate heading differences ----
  dtheta <- array(
    rnorm(m * m * T * L, sd = 0.5),
    dim = c(m, m, T, L)
  )
  
  ## ---- Simulate counts ----
  Y <- matrix(0, m, T)
  
  for (b in (L + 1):T) {
    for (i in 1:m) {
      
      lambda <- mu_true[i]
      
      for (j in 1:m) {
        if (j != i) {
          for (ell in 1:L) {
            lambda <- lambda +
              alpha_true[i, j] *
              Y[j, b - ell] *
              exp(-beta_true * ell) *
              exp(-(dtheta[i, j, b, ell]^2) / (2 * sigma_true^2))
          }
        }
      }
      
      Y[i, b] <- rpois(1, lambda)
    }
  }
  
  ## ---- Optimize ----
  fit <- tryCatch(
    optim(
      par     = init_par,
      fn      = neg_loglik,
      Y       = Y,
      dtheta  = dtheta,
      m       = m,
      L       = L,
      method  = "BFGS",
      hessian = TRUE,
      control = list(maxit = 500)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit) || fit$convergence != 0) {
    converged[sim] <- FALSE
    next
  }
  
  ## ---- Variance-covariance ----
  vcov_mat <- tryCatch(solve(fit$hessian), error = function(e) NULL)
  if (is.null(vcov_mat)) {
    converged[sim] <- FALSE
    next
  }
  
  converged[sim] <- TRUE
  se_par <- sqrt(diag(vcov_mat))
  
  ## ---- Beta CI ----
  beta_hat <- exp(fit$par[m + n_alpha + 1])
  beta_se  <- se_par[m + n_alpha + 1]
  beta_ci  <- exp(log(beta_hat) + c(-1, 1) * z_alpha * beta_se)
  covered[sim, "beta"] <- (beta_ci[1] <= beta_true) & (beta_true <= beta_ci[2])
  
  ## ---- Sigma CI ----
  sigma_hat <- exp(fit$par[m + n_alpha + 2])
  sigma_se  <- se_par[m + n_alpha + 2]
  sigma_ci  <- exp(log(sigma_hat) + c(-1, 1) * z_alpha * sigma_se)
  covered[sim, "sigma"] <- (sigma_ci[1] <= sigma_true) & (sigma_true <= sigma_ci[2])
  
  ## ---- Alpha CIs for selected pairs ----
  target_pairs <- list(c(2,1), c(3,2), c(4,3), c(5,4))
  
  for (pair in target_pairs) {
    i <- pair[1]; j <- pair[2]
    pname <- paste0("alpha[", i, ",", j, "]")
    k     <- idx_map[i, j]
    
    eta_k    <- fit$par[m + k]
    se_eta_k <- se_par[m + k]
    
    ci_lo <- exp(eta_k - z_alpha * se_eta_k)
    ci_hi <- exp(eta_k + z_alpha * se_eta_k)
    
    covered[sim, pname] <- (ci_lo <= alpha_true[i, j]) & (alpha_true[i, j] <= ci_hi)
  }
  print(sim)
}

#############################################
# SUMMARISE RESULTS
#############################################

n_converged <- sum(converged)
cat("\n=== MONTE CARLO COVERAGE RESULTS ===\n")
cat("Total simulations :", N_sim, "\n")
cat("Converged         :", n_converged, "\n\n")

coverage_rates <- colMeans(covered, na.rm = TRUE)

results_df <- data.frame(
  Parameter    = param_names,
  True_Value   = round(true_vals, 3),
  Coverage     = round(coverage_rates, 3),
  Coverage_Pct = paste0(round(coverage_rates * 100, 1), "%")
)

print(results_df, row.names = FALSE)

cat("\nNominal coverage level: 95%\n")