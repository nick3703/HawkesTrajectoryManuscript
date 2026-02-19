# Marked Multivariate Hawkes for Coordinated Movement - MULTIPLE LEADERS VERSION
# -------------------------------------------------------------------------------
# Modified to test detection of multiple leaders in different scenarios:
#   Scenario 1: Single dominant leader (baseline)
#   Scenario 2: Two co-leaders with equal influence
#   Scenario 3: Hierarchical leadership (primary + secondary leader)
#   Scenario 4: Dynamic leadership (leaders influence different subgroups)
# -------------------------------------------------------------------------------

library(circular)
library(ggplot2)
library(dplyr)

# ---------------------------
# 0) Enhanced Simulation with Multiple Leader Scenarios
# ---------------------------

simulate_pack_multi_leader <- function(n_agents, n_steps, dt, 
                                       scenario = c("single", "co_leaders", 
                                                    "hierarchical", "subgroups"),
                                       leader_align_prob = 0.7,
                                       speed_mean = 1.0, speed_sd = 0.3) {
  
  scenario <- match.arg(scenario)
  
  x <- matrix(0, nrow = n_steps, ncol = n_agents)
  y <- matrix(0, nrow = n_steps, ncol = n_agents)
  theta <- matrix(NA_real_, nrow = n_steps, ncol = n_agents)
  v <- matrix(NA_real_, nrow = n_steps, ncol = n_agents)
  
  # Initial headings and speeds
  theta[1, ] <- runif(n_agents, -pi, pi)
  v[1, ] <- pmax(0.05, rnorm(n_agents, speed_mean, speed_sd))
  
  for (k in 2:n_steps) {
    
    if (scenario == "single") {
      # ========================================
      # Scenario 1: Single leader (agent 1)
      # ========================================
      # Leader evolves independently
      theta[k, 1] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 1]), kappa = 4)
      v[k, 1] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # All followers align to agent 1
      for (j in 2:n_agents) {
        if (runif(1) < leader_align_prob) {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k, 1]), kappa = 3)
        } else {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, j]), kappa = 3)
        }
        v[k, j] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      }
      
    } else if (scenario == "co_leaders") {
      # ========================================
      # Scenario 2: Two co-leaders (agents 1 & 2)
      # ========================================
      # Both leaders evolve somewhat independently but with slight coordination
      theta[k, 1] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 1]), kappa = 4)
      v[k, 1] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # Leader 2 has some tendency to align with leader 1 (weak coupling)
      if (runif(1) < 0.3) {
        theta[k, 2] <- circular::rvonmises(1, mu = circular::circular(theta[k, 1]), kappa = 2)
      } else {
        theta[k, 2] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 2]), kappa = 4)
      }
      v[k, 2] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # Followers (3-5) align to EITHER leader 1 or 2 with equal probability
      for (j in 3:n_agents) {
        if (runif(1) < leader_align_prob) {
          # Choose which leader to follow
          leader_to_follow <- sample(c(1, 2), 1)
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k, leader_to_follow]), kappa = 3)
        } else {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, j]), kappa = 3)
        }
        v[k, j] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      }
      
    } else if (scenario == "hierarchical") {
      # ========================================
      # Scenario 3: Hierarchical (1 is primary, 2 is secondary)
      # ========================================
      # Primary leader (agent 1)
      theta[k, 1] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 1]), kappa = 4)
      v[k, 1] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # Secondary leader (agent 2) follows agent 1 with high probability
      if (runif(1) < 0.8) {
        theta[k, 2] <- circular::rvonmises(1, mu = circular::circular(theta[k, 1]), kappa = 3)
      } else {
        theta[k, 2] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 2]), kappa = 3)
      }
      v[k, 2] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # Followers 3-4 align primarily to agent 2 (secondary leader)
      for (j in 3:4) {
        if (runif(1) < leader_align_prob) {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k, 2]), kappa = 3)
        } else {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, j]), kappa = 3)
        }
        v[k, j] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      }
      
      # Follower 5 aligns to agent 1 directly
      if (runif(1) < leader_align_prob) {
        theta[k, 5] <- circular::rvonmises(1, mu = circular::circular(theta[k, 1]), kappa = 3)
      } else {
        theta[k, 5] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 5]), kappa = 3)
      }
      v[k, 5] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
    } else if (scenario == "subgroups") {
      # ========================================
      # Scenario 4: Subgroup leaders (1 leads 2-3, 4 leads 5)
      # ========================================
      # Leader 1 (for subgroup 2-3)
      theta[k, 1] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 1]), kappa = 4)
      v[k, 1] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # Followers 2-3 align to agent 1
      for (j in 2:3) {
        if (runif(1) < leader_align_prob) {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k, 1]), kappa = 3)
        } else {
          theta[k, j] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, j]), kappa = 3)
        }
        v[k, j] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      }
      
      # Leader 4 (for subgroup containing just 5)
      theta[k, 4] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 4]), kappa = 4)
      v[k, 4] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
      
      # Follower 5 aligns to agent 4
      if (runif(1) < leader_align_prob) {
        theta[k, 5] <- circular::rvonmises(1, mu = circular::circular(theta[k, 4]), kappa = 3)
      } else {
        theta[k, 5] <- circular::rvonmises(1, mu = circular::circular(theta[k-1, 5]), kappa = 3)
      }
      v[k, 5] <- pmax(0.05, rnorm(1, speed_mean, speed_sd))
    }
    
    # Update positions for all agents
    x[k, ] <- x[k-1, ] + as.numeric(v[k, ] * cos(theta[k, ]) * dt)
    y[k, ] <- y[k-1, ] + as.numeric(v[k, ] * sin(theta[k, ]) * dt)
  }
  
  # Assemble dataframe
  t <- seq(0, by = dt, length.out = n_steps)
  rows <- do.call(rbind, lapply(1:n_agents, function(j) {
    cbind(id = j, t = t, x = x[, j], y = y[, j])
  }))
  df <- as.data.frame(rows)
  df$id <- factor(df$id, levels = 1:n_agents)
  df
}

# ---------------------------
# Keep all your existing helper functions unchanged
# ---------------------------

# Helper: minimal angular difference in (-pi, pi]
ang_diff <- function(a, b) {
  d <- a - b
  d <- (d + pi) %% (2*pi) - pi
  d
}

eventize <- function(df, dt = 10, v_min = 0.2) {
  split_df <- split(df, df$id)
  ev_list <- lapply(split_df, function(d) {
    d <- d[order(d$t), ]
    dx <- c(NA, diff(d$x))
    dy <- c(NA, diff(d$y))
    v <- sqrt(dx^2 + dy^2) / dt
    th <- atan2(dy, dx)
    keep <- which(!is.na(v) & v >= v_min)
    if (length(keep) == 0) return(NULL)
    data.frame(id = d$id[keep], t = d$t[keep], x = d$x[keep], y = d$y[keep],
               v = v[keep], theta = th[keep])
  })
  do.call(rbind, ev_list)
}

build_pairs <- function(ev, T_max) {
  n <- nrow(ev)
  starts <- integer(n)
  j <- 1
  for (i in 1:n) {
    ti <- ev$t[i]
    while (j < i && ti - ev$t[j] > T_max) j <- j + 1
    starts[i] <- j
  }
  i_vec <- integer()
  k_vec <- integer()
  if (n > 1) {
    for (i in 2:n) {
      if (starts[i] <= i-1) {
        ks <- starts[i]:(i-1)
        i_vec <- c(i_vec, rep(i, length(ks)))
        k_vec <- c(k_vec, ks)
      }
    }
  }
  if (length(i_vec) == 0) return(data.frame())
  dt <- ev$t[i_vec] - ev$t[k_vec]
  dx <- ev$x[i_vec] - ev$x[k_vec]
  dy <- ev$y[i_vec] - ev$y[k_vec]
  dist <- sqrt(dx^2 + dy^2)
  dv <- ev$v[i_vec] - ev$v[k_vec]
  dth <- ang_diff(ev$theta[i_vec], ev$theta[k_vec])
  data.frame(i_eid = ev$eid[i_vec], k_eid = ev$eid[k_vec],
             id_i = ev$id[i_vec], id_k = ev$id[k_vec],
             dt = dt, dist = dist, dv = dv, dtheta = dth)
}

# [Keep all your kernel functions and estimation code exactly as is]
# [I'm including the key estimation wrapper here]

fit_hawkes <- function(events, pairs, S, agents, T_obs, T_max, beta, T_end) {
  
  m <- length(agents)
  I_map <- setNames(seq_along(agents), agents)
  
  # Parameter setup
  mu_idx <- 1:m
  alpha_map <- which(row(matrix(1, m, m)) != col(matrix(1, m, m)), arr.ind = TRUE)
  n_alpha <- nrow(alpha_map)
  alpha_idx <- (m + 1):(m + n_alpha)
  
  build_Alpha <- function(theta) {
    Alpha <- matrix(0, nrow = m, ncol = m)
    Alpha[alpha_map] <- theta[alpha_idx]
    Alpha
  }
  
  lambda_at_events <- function(theta) {
    mu <- theta[mu_idx]
    Alpha <- build_Alpha(theta)
    id_i_row <- I_map[as.character(events$id)]
    lam <- mu[id_i_row]
    for (j in 1:m) lam <- lam + Alpha[id_i_row, j] * S[, j]
    pmax(lam, .Machine$double.eps)
  }
  
  by_src <- split(events$t, events$id)
  G_fun <- function(u) {
    u_cap <- pmin(u, T_max)
    1 - exp(-beta * u_cap)
  }
  G_j <- sapply(agents, function(a) sum(G_fun(T_end - by_src[[a]])))
  
  loglik <- function(theta) {
    theta <- pmax(theta, 1e-10)
    lam <- lambda_at_events(theta)
    ll1 <- sum(log(lam))
    mu <- theta[mu_idx]
    Alpha <- build_Alpha(theta)
    G_mat <- matrix(rep(G_j, each = m), nrow = m, byrow = TRUE)
    comp <- sum(mu) * T_obs + sum(Alpha * G_mat)
    ll1 - comp
  }
  
  # Initialize
  mu0 <- rep(nrow(events) / (m * T_obs + 1e-6), m)
  alpha0 <- rep(0.005, n_alpha)
  theta0 <- c(mu0, alpha0)
  lower <- rep(0, length(theta0))
  upper <- rep(Inf, length(theta0))
  
  opt <- optim(theta0, fn = function(th) -loglik(th), method = "L-BFGS-B",
               lower = lower, upper = upper, control = list(maxit = 200))
  
  theta_hat <- pmax(opt$par, 0)
  mu_hat <- theta_hat[mu_idx]
  Alpha_hat <- matrix(0, nrow = m, ncol = m)
  Alpha_hat[alpha_map] <- theta_hat[alpha_idx]
  rownames(Alpha_hat) <- paste0("target_", agents)
  colnames(Alpha_hat) <- paste0("source_", agents)
  
  list(mu = mu_hat, Alpha = Alpha_hat, theta = theta_hat, opt = opt)
}

# ---------------------------
# Main Simulation Loop for Multiple Scenarios
# ---------------------------

run_scenario_simulation <- function(scenario_name, n_reps = 500, 
                                    leader_align_prob = 0.7) {
  
  cat(sprintf("\n========== Running Scenario: %s ==========\n", scenario_name))
  
  n_agents <- 5
  T_total <- 30 * 60
  dt <- 10
  n_steps <- T_total / dt
  T_max <- 120
  
  # Kernel parameters
  beta <- 1/30
  sigma_d <- 25
  sigma_v <- 0.5
  kappa <- 4
  
  phi_t <- function(dt, beta) ifelse(dt > 0 & dt <= T_max, beta * exp(-beta * dt), 0)
  phi_d <- function(r, sd) exp(-(r^2)/(2 * sd^2))
  phi_v <- function(dv, sv) exp(-(dv^2)/(2 * sv^2))
  phi_th <- function(dth, kappa) exp(kappa * cos(dth))
  
  results_net <- matrix(NA, nrow = n_reps, ncol = n_agents)
  results_alpha <- array(NA, dim = c(n_agents, n_agents, n_reps))
  
  for (rep in 1:n_reps) {
    if (rep %% 100 == 0) cat(sprintf("  Rep %d/%d\n", rep, n_reps))
    
    # Simulate
    df <- simulate_pack_multi_leader(n_agents, n_steps, dt, 
                                     scenario = scenario_name,
                                     leader_align_prob = leader_align_prob)
    
    # Eventize
    events <- eventize(df, dt = dt, v_min = 0.2)
    if (nrow(events) < 10) next  # Skip if too few events
    
    rownames(events) <- NULL
    ord <- order(events$t)
    events <- events[ord, ]
    events$eid <- seq_len(nrow(events))
    
    # Build pairs
    pairs <- build_pairs(events, T_max)
    if (nrow(pairs) == 0) next
    
    # Compute weights
    pairs$w <- with(pairs, phi_t(dt, beta) * phi_d(dist, sigma_d) * 
                      phi_v(dv, sigma_v) * phi_th(dtheta, kappa))
    
    # Build S matrix
    agents <- levels(events$id)
    I_map <- setNames(seq_along(agents), agents)
    S <- matrix(0, nrow = nrow(events), ncol = length(agents))
    colnames(S) <- agents
    
    key <- paste(pairs$i_eid, pairs$id_k, sep = "|")
    agg <- tapply(pairs$w, key, sum)
    keys <- strsplit(names(agg), "\\|")
    i_idx <- as.integer(sapply(keys, `[[`, 1))
    idk <- sapply(keys, `[[`, 2)
    col_idx <- I_map[as.character(idk)]
    S[cbind(i_idx, col_idx)] <- as.numeric(agg)
    
    # Fit model
    T_start <- min(events$t)
    T_end <- max(events$t)
    T_obs <- T_end - T_start
    
    fit <- fit_hawkes(events, pairs, S, agents, T_obs, T_max, beta, T_end)
    
    # Store results
    L_out <- colSums(fit$Alpha)
    F_in <- rowSums(fit$Alpha)
    L_net <- L_out - F_in
    
    results_net[rep, ] <- L_net
    results_alpha[,, rep] <- fit$Alpha
  }
  
  # Analyze results
  valid_rows <- complete.cases(results_net)
  results_net <- results_net[valid_rows, ]
  
  # Who is detected as leader most often?
  max_leader <- apply(results_net, 1, which.max)
  detection_freq <- table(max_leader) / sum(valid_rows)
  
  # Mean alpha matrix
  mean_alpha <- apply(results_alpha[,, valid_rows], c(1,2), mean, na.rm = TRUE)
  rownames(mean_alpha) <- paste0("Target_", 1:n_agents)
  colnames(mean_alpha) <- paste0("Source_", 1:n_agents)
  
  # Mean net leadership
  mean_net <- colMeans(results_net, na.rm = TRUE)
  
  list(
    scenario = scenario_name,
    detection_freq = detection_freq,
    mean_alpha = mean_alpha,
    mean_net = mean_net,
    results_net = results_net,
    n_valid = sum(valid_rows)
  )
}

# ---------------------------
# Run all scenarios
# ---------------------------

set.seed(123)

results_single <- run_scenario_simulation("single", n_reps = 500, leader_align_prob = 0.7)
results_co <- run_scenario_simulation("co_leaders", n_reps = 500, leader_align_prob = 0.7)
results_hier <- run_scenario_simulation("hierarchical", n_reps = 500, leader_align_prob = 0.7)
results_sub <- run_scenario_simulation("subgroups", n_reps = 500, leader_align_prob = 0.7)

# ---------------------------
# Summary outputs
# ---------------------------

cat("\n\n========== SUMMARY OF ALL SCENARIOS ==========\n\n")

print_scenario_summary <- function(res) {
  cat(sprintf("\nScenario: %s (n=%d valid reps)\n", res$scenario, res$n_valid))
  cat("\nLeader Detection Frequency:\n")
  print(round(res$detection_freq, 3))
  cat("\nMean Net Leadership:\n")
  print(round(res$mean_net, 4))
  cat("\nMean Alpha Matrix:\n")
  print(round(res$mean_alpha, 4))
  cat("\n" , paste(rep("-", 60), collapse = ""), "\n")
}

print_scenario_summary(results_single)
print_scenario_summary(results_co)
print_scenario_summary(results_hier)
print_scenario_summary(results_sub)

# ---------------------------
# Visualization
# ---------------------------

# Plot net leadership distributions by scenario
plot_data <- rbind(
  data.frame(scenario = "Single Leader", 
             entity = rep(1:5, each = nrow(results_single$results_net)),
             net_leadership = c(results_single$results_net)),
  data.frame(scenario = "Co-Leaders", 
             entity = rep(1:5, each = nrow(results_co$results_net)),
             net_leadership = c(results_co$results_net)),
  data.frame(scenario = "Hierarchical", 
             entity = rep(1:5, each = nrow(results_hier$results_net)),
             net_leadership = c(results_hier$results_net)),
  data.frame(scenario = "Subgroups", 
             entity = rep(1:5, each = nrow(results_sub$results_net)),
             net_leadership = c(results_sub$results_net))
)

ggplot(plot_data, aes(x = factor(entity), y = net_leadership, fill = factor(entity))) +
  geom_boxplot() +
  facet_wrap(~scenario, scales = "free_y") +
  theme_bw() +
  labs(x = "Entity ID", y = "Net Leadership Score", fill = "Entity",
       title = "Net Leadership Distribution Across Scenarios") +
  theme(legend.position = "bottom")

# Heatmap comparison of mean alpha matrices
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))

plot_alpha_heatmap <- function(alpha_mat, title) {
  image(1:ncol(alpha_mat), 1:nrow(alpha_mat), t(alpha_mat[nrow(alpha_mat):1,]),
        col = hcl.colors(20, "YlOrRd", rev = FALSE),
        xlab = "Source", ylab = "Target", main = title,
        axes = FALSE)
  axis(1, at = 1:ncol(alpha_mat), labels = 1:ncol(alpha_mat))
  axis(2, at = 1:nrow(alpha_mat), labels = nrow(alpha_mat):1)
  box()
}

plot_alpha_heatmap(results_single$mean_alpha, "Single Leader")
plot_alpha_heatmap(results_co$mean_alpha, "Co-Leaders")
plot_alpha_heatmap(results_hier$mean_alpha, "Hierarchical")
plot_alpha_heatmap(results_sub$mean_alpha, "Subgroups")