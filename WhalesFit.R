#############################################
# DISCRETIZED HAWKES TRAJECTORY
# Whale Movement Example
#############################################

rm(list = ls())
set.seed(123)

#############################################
# 1. LIBRARIES
#############################################
library(data.table)
library(Matrix)
library(glmnet)
library(geosphere)

#############################################
# 2. LOAD AND PREPROCESS DATA
#############################################

dat <- fread("Whales-data.csv")

# Create timestamp
dat[, time := as.POSIXct(
  sprintf("20%02d-%02d-%02d %02d:%02d:%02d",
          Year, Month, Day, Hour, Min, Sec),
  tz = "UTC"
)]

# Order data
setorder(dat, ID, time)


dat<-dat %>% group_by(ID) %>%
  filter(max(abs(X))<125 & min(abs(X))>118)%>%
  filter(max(Y)<37&min(Y)>32)%>%
  ungroup()
dat <- as.data.table(dat)

# Encode whales
dat[, whale := as.integer(factor(ID))]
m <- length(unique(dat$whale))





#############################################
# 3. COMPUTE SPEED AND HEADING
#############################################
dat[, `:=`(
  lon = X,
  lat = Y
)]

dat[, `:=`(
  lon_prev = shift(lon),
  lat_prev = shift(lat),
  time_prev = shift(time)
), by = whale]

dat[, dt := as.numeric(difftime(time, time_prev, units = "secs"))]

dat[, dist := distHaversine(
  cbind(lon_prev, lat_prev),
  cbind(lon, lat)
)]

dat[, speed := dist / dt]

dat[, heading := bearing(
  cbind(lon_prev, lat_prev),
  cbind(lon, lat)
)]

dat <- dat[!is.na(speed) & dt > 0]

#############################################
# 4. DISCRETIZE TIME
#############################################
Delta <-  3600   # 1 hour bins
t0 <- min(dat$time)
dat[, bin := as.integer(floor(as.numeric(difftime(time, t0, units="secs")) / Delta))]

B <- max(dat$bin) + 1

#############################################
# 5. CREATE COUNT MATRIX Y_{i,b}
#############################################
Y <- sparseMatrix(
  i = dat$whale,
  j = dat$bin+1 ,
  x=1,
  dims = c(m, B)
)

#############################################
# 6. KERNEL SETUP (TEMPORAL + HEADING)
#############################################
L <- 10                       # number of lags
beta <- 0.8                   # temporal decay
sigma_theta <- 0.5            # angular bandwidth

time_kernel <- exp(-beta * (1:L))
time_kernel <- time_kernel / sum(time_kernel)
#############################################
# 7. HEADING MATRIX (VECTOR FORM)
#############################################
heading_mat <- matrix(NA, m, B)

for (i in 1:m) {
  idx <- which(dat$whale == i)
  heading_mat[i, dat$bin[idx] + 1] <- dat$heading[idx]
}

#############################################
# 8. BUILD SPARSE DESIGN MATRICES (SAFE)
#############################################
X_list <- vector("list", m)
y_list <- vector("list", m)

for (i in 1:m) {
  
  y_i <- as.numeric(Y[i, ])
  
  X_cols <- list()
  col_idx <- 1
  
  for (j in 1:m) {
    if (j == i) next
    
    x_ij <- numeric(B)
    y_conv <- numeric(B)
    for (l in 1:L) {
      y_conv[(l+1):B] <- y_conv[(l+1):B] +
        exp(-beta * l) * Y[j, 1:(B-l)]
    }
    
    # heading weight (no lag dependence)
    theta_diff <- heading_mat[i, ] - heading_mat[j, ]
    theta_weight <- exp(-(theta_diff^2) / (2 * sigma_theta^2))
    theta_weight[is.na(theta_weight)] <- 1   # IMPORTANT
    
    x_ij <- y_conv * theta_weight
    X_cols[[col_idx]] <- x_ij
    col_idx <- col_idx + 1
  }
  
  X_i <- do.call(cbind, X_cols)
  
  X_list[[i]] <- Matrix(X_i, sparse = TRUE)
  y_list[[i]] <- y_i
  if (i %% 10 == 0) cat("Finished target", i, "\n")
}

zero_matrices <- sapply(X_list, function(mat) all(max(mat) > .0005))

col_norms <- sqrt(colSums(X_i^2))
keep_cols <- col_norms > 1e-6

if (!any(keep_cols)) next

X_i <- X_i[, keep_cols, drop = FALSE]


fit_hawkes_unpenalized <- function(beta, sigma_theta,
                                   Y, heading_mat, m, B, L) {
  
  total_loglik <- 0
  mu_hat <- numeric(m)
  
  for (i in 1:m) {
    
    y_i <- as.numeric(Y[i, ])
    X_cols <- list()
    
    for (j in 1:m) {
      if (j == i) next
      
      x_ij <- numeric(B)
      
      for (l in 1:L) {
        y_lag <- numeric(B)
        y_lag[(l + 1):B] <- Y[j, 1:(B - l)]
        
        theta_i <- heading_mat[i, ]
        theta_j <- heading_mat[j, ]
        
        theta_diff <- theta_i - c(rep(NA, l), theta_j[1:(B - l)])
        theta_weight <- exp(-(theta_diff^2) / (2 * sigma_theta^2))
        theta_weight[is.na(theta_weight)] <- 0
        
        x_ij <- x_ij +
          y_lag *
          exp(-beta * l) *
          theta_weight
      }
      
      X_cols[[length(X_cols) + 1]] <- x_ij
    }
    
    X_i <- do.call(cbind, X_cols)
    if (all(X_i == 0)) next
    
    fit <- glm(y_i~X_i,
      family = "poisson",
    )
    
    eta <- as.numeric(fitted(fit, type = "link"))
    total_loglik <- total_loglik +
      sum(y_i * eta - exp(eta))
    
    mu_hat[i] <- coef(fit)[1]
  }
  
  total_loglik
}


objective_kernel <- function(par, Y, heading_mat, m, B, L) {
  beta <- 0
  sigma_theta <- exp(par[1])
  
  ll <- fit_hawkes_unpenalized(
    beta = beta,
    sigma_theta = sigma_theta,
    Y = Y,
    heading_mat = heading_mat,
    m = m,
    B = B,
    L = L
  )
  
  -ll  # optim minimizes
}

# Initial values
par_init <- log(c(sigma_theta = 0.5))

opt <- optim(
  par = par_init,
  fn = objective_kernel,
  Y = Y,
  heading_mat = heading_mat,
  m = m,
  B = B,
  L = L,
  method = "BFGS",
  control = list(maxit = 50)
)

#beta_hat  <- exp(opt$par[1])
sigma_hat <- exp(opt$par[1])

#cat("Estimated beta:", beta_hat, "\n")
cat("Estimated sigma_theta:", sigma_hat, "\n")



fit_hawkes_given_kernels <- function(beta, sigma_theta,
                                     Y, heading_mat, m, B, L) {
  
  alpha_hat <- matrix(0, m, m)
  mu_hat <- numeric(m)
  total_loglik <- 0
  
  for (i in 1:m) {
    
    y_i <- as.numeric(Y[i, ])
    X_cols <- list()
    col_idx <- 1
    
    for (j in 1:m) {
      if (j == i) next
      
      x_ij <- numeric(B)
      
      for (l in 1:L) {
        y_lag <- numeric(B)
        y_lag[(l+1):B] <- Y[j, 1:(B-l)]
        
        theta_i <- heading_mat[i, ]
        theta_j <- heading_mat[j, ]
        
        theta_diff <- theta_i - c(rep(NA, l), theta_j[1:(B-l)])
        theta_weight <- exp(-(theta_diff^2) / (2 * sigma_theta^2))
        theta_weight[is.na(theta_weight)] <- 0
        
        x_ij <- x_ij + y_lag * exp(-beta * l) * theta_weight
      }
      
      X_cols[[col_idx]] <- x_ij
      col_idx <- col_idx + 1
    }
    
    X_i <- Matrix(do.call(cbind, X_cols), sparse = TRUE)
    
    if (all(X_i == 0)) next
    
    fit <- glmnet(
      x = X_i,
      y = y_i,
      family = "poisson",
      alpha = 1,
      lambda = 1e-6,   # FIXED penalty during kernel optimization
      standardize = FALSE
    )
    
    coef_i <- as.numeric(coef(fit))
    mu_hat[i] <- coef_i[1]
    
    idx <- 2
    for (j in 1:m) {
      if (j == i) next
      alpha_hat[i, j] <- coef_i[idx]
      idx <- idx + 1
    }
    
    # contribution to log-likelihood
    eta <- as.numeric(X_i %*% coef_i[-1] + coef_i[1])
    total_loglik <- total_loglik +
      sum(y_i * eta - exp(eta))
  }
  
  list(
    loglik = total_loglik,
    alpha = alpha_hat,
    mu = mu_hat
  )
}


objective <- function(par) {
  beta <- exp(par[1])
  sigma_theta <- exp(par[2])
  
  fit <- fit_hawkes_given_kernels(
    beta = beta,
    sigma_theta = sigma_theta,
    Y = Y,
    heading_mat = heading_mat,
    m = m,
    B = B,
    L = L
  )
  
  # Negative because optim minimizes
  -fit$loglik
}

# Initial guesses
par0 <- log(c(beta = 0.3, sigma_theta = 0.2))

opt <- optim(
  par = par0,
  fn = objective,
  method = "Nelder-Mead",
  control = list(maxit = 50)
)

beta_hat <- exp(opt$par[1])
sigma_theta_hat <- exp(opt$par[2])

cat("Estimated beta:", beta_hat, "\n")
cat("Estimated sigma_theta:", sigma_theta_hat, "\n")





#############################################
# 9. FIT SPARSE POISSON GLM (ROW-BY-ROW)
#############################################
alpha_hat <- matrix(0, m, m)
mu_hat <- numeric(m)

for (i in 1:m) {
  
  X_i <- X_list[[i]]
  y_i <- y_list[[i]]
  
  fit <- cv.glmnet(
    x = X_i,
    y = y_i,
    family = "poisson",
    alpha = 1,            # LASSO for sparsity
    standardize = FALSE
  )
  
  beta_hat <- as.numeric(coef(fit, s = "lambda.min"))
  mu_hat[i] <- beta_hat[1]
  
  idx <- 2
  for (j in 1:m) {
    if (j == i) next
    alpha_hat[i, j] <- beta_hat[idx]
    idx <- idx + 1
  }
  print(i)
}

#############################################
# 10. OUTPUT RESULTS
#############################################
cat("\nEstimated alpha matrix:\n")
print(round(alpha_hat, 3))

cat("\nEstimated background rates:\n")
print(round(mu_hat, 4))
