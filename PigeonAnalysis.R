
rm(list = ls()); set.seed(123)
library(data.table)
library(sf)

## ---------------------------------------------------------------------
## 0. CONFIG
## ---------------------------------------------------------------------
CFG <- list(
  data_dir   = "Data",                 ## VERIFY: folder holding all the CSVs
  pattern    = "^01_.*\\.csv$",     ## Dryad naming: 01_<flock>_<bird>_<release>.csv
  utm_epsg   = 32630,               ## UTM 30N -- Oxford/Wytham, UK (~1.3 W). VERIFY
  
  speed_kmh_fly = 15,               ## SPEED above this (km/h) = in flight. TUNE
  min_fly_sec   = 30,               ## require a flight segment at least this long
  
  bin_sec    = 0.4,                 ## 5 Hz data; bins of ~0.4 s. KEY KNOB (0.2-1.0). TUNE
  L          = 8,                   ## lags; L*bin_sec ~ 3 s memory (delays ~0.3 s). TUNE
  beta       = 1.5,                 ## temporal decay over sub-second lags. TUNE
  kappa      = 6,                   ## von Mises heading concentration (tighter than baboons)
  sigma_s    = 15,                  ## spatial bandwidth, METRES (flock spans ~tens of m). TUNE
  d_max      = 60,                  ## proximity gate, METRES. TUNE
  min_events = 20                   ## min movement events for a bird to be a target
)

## =====================================================================
## 1. MERGE PIPELINE  -- read every CSV, recover identity from FILENAME
## =====================================================================
files <- list.files(CFG$data_dir, pattern = CFG$pattern, full.names = TRUE,
                    recursive = TRUE)
stopifnot(length(files) > 0)
cat("found", length(files), "CSV files\n")

parse_one <- function(path) {
  fn   <- sub("\\.csv$", "", basename(path))
  part <- strsplit(fn, "_")[[1]]            ## c("01", flock, bird, release)
  flock <- part[2]; bird <- part[3]; release <- part[4]
  d <- fread(path)
  setnames(d, trimws(names(d)))             ## strip stray spaces in headers
  d <- d[VALID == "FIXED"]                   ## keep good fixes only
  if (nrow(d) == 0) return(NULL)
  ## timestamp = UTC DATE + UTC TIME + MS
  d[, time := as.POSIXct(paste(`UTC DATE`, `UTC TIME`), tz = "UTC",
                         format = "%Y/%m/%d %H:%M:%S") + MS/1000]
  ## signed coordinates (W -> negative lon, S -> negative lat)
  d[, lat := fifelse(`N/S` == "S", -1, 1) * LATITUDE]
  d[, lon := fifelse(`E/W` == "W", -1, 1) * LONGITUDE]
  d[, .(flock, bird, release, time, lon, lat,
        speed = SPEED, gps_heading = HEADING)]
}
all_raw <- rbindlist(lapply(files, parse_one), use.names = TRUE)
setorder(all_raw, flock, release, bird, time)
all_raw[, frkey := paste(flock, release, sep = "_")]   ## one flock-flight = analysis unit
cat("merged:", nrow(all_raw), "fixes across",
    uniqueN(all_raw$frkey), "flock-releases,",
    uniqueN(all_raw[, .(flock, bird)]), "distinct birds\n")


## ---------- 2. PROJECT to metres (fast matrix path; no geometry objects) ----------
xy <- sf::sf_project(from = "EPSG:4326",
                     to   = paste0("EPSG:", CFG$utm_epsg),
                     pts  = as.matrix(all_raw[, .(lon, lat)]))
all_raw[, `:=`(x = xy[, 1], y = xy[, 2])]

## =====================================================================
##  ESTIMATOR (identity-link Poisson; L1+BIC+relaxed refit; Fisher SE)
##  -- unchanged from the validated whale/baboon estimator.
## =====================================================================
## E1
fit_target_idlink <- function(S, y, P, eps = 1e-8, theta0 = NULL) {
  p <- ncol(S); A <- cbind(1, S)
  if (is.null(theta0)) theta0 <- c(max(mean(y), eps), rep(0, p))
  negll <- function(th){lam <- as.numeric(A %*% th); -sum(y*log(pmax(lam,1e-12))-lam)+P*sum(th[-1])}
  grad  <- function(th){lam <- pmax(as.numeric(A %*% th),1e-12)
  g <- -as.numeric(crossprod(A, y/lam-1)); g[-1] <- g[-1]+P; g}
  opt <- optim(theta0, negll, grad, method = "L-BFGS-B",
               lower = c(eps, rep(0, p)), control = list(maxit = 1000, factr = 1e7))
  lam <- pmax(as.numeric(A %*% opt$par), 1e-12)
  list(mu = opt$par[1], alpha = opt$par[-1], theta = opt$par,
       loglik = sum(y*log(lam)-lam))
}
## E2
fit_target_path <- function(S, y, ngrid = 20, tol = 1e-6, relax = TRUE) {
  n <- length(y); mu0 <- max(mean(y), 1e-8)
  Pmax  <- max(c(as.numeric(crossprod(S, y/mu0 - 1)), 1))
  Pgrid <- Pmax * 10^seq(0, -3, length.out = ngrid)
  best <- NULL; bestic <- Inf; theta0 <- NULL
  for (P in Pgrid) {
    f <- fit_target_idlink(S, y, P, theta0 = theta0); theta0 <- f$theta
    ic <- -2*f$loglik + log(n)*(1 + sum(f$alpha > tol))
    if (ic < bestic) { bestic <- ic; best <- f; best$P <- P }
  }
  if (relax) {
    keep <- which(best$alpha > tol)
    fr <- fit_target_idlink(S[, keep, drop = FALSE], y, P = 0)
    a <- numeric(ncol(S)); if (length(keep)) a[keep] <- fr$alpha
    best$mu <- fr$mu; best$alpha <- a
  }
  best
}
## E3
target_cov <- function(S, y, mu, alpha, tol = 1e-6) {
  A <- cbind(1, S); act <- c(TRUE, alpha > tol); Aa <- A[, act, drop = FALSE]
  lam <- pmax(as.numeric(Aa %*% c(mu, alpha)[act]), 1e-12)
  info <- crossprod(Aa, Aa / lam)
  V <- tryCatch(solve(info), error = function(e) solve(info + 1e-10*diag(ncol(info))))
  V[-1, -1, drop = FALSE]
}

## =====================================================================
## 3. PER FLOCK-RELEASE ANALYSIS
##    Returns the alpha matrix, leadership, and directional validation
##    for one flight. Pooling across flights happens in section 4.
## =====================================================================
ang_diff <- function(a, b) abs(atan2(sin(a - b), cos(a - b)))

analyze_flight <- function(D) {
  ## D: rows for ONE flock-release (multiple birds). Isolate the flight,
  ## bin, build events+marks, fit, validate.
  birds <- sort(unique(D$bird))
  if (length(birds) < 2) return(NULL)
  
  ## --- isolate the flight segment: keep the contiguous window where the
  ##     flock is actually moving (>= speed threshold for >= min_fly_sec) ---
  D[, fly := speed >= CFG$speed_kmh_fly]
  flo <- D[fly == TRUE]
  if (nrow(flo) == 0) return(NULL)
  tlo <- min(flo$time); thi <- max(flo$time)
  if (as.numeric(thi - tlo, units = "secs") < CFG$min_fly_sec) return(NULL)
  D <- D[time >= tlo & time <= thi]
  
  ## --- common sub-second grid ---
  t0 <- as.numeric(min(D$time))
  D[, bin := as.integer(floor((as.numeric(time) - t0) / CFG$bin_sec))]
  reg <- D[, .(x = mean(x), y = mean(y)), by = .(bird, bin)]
  setorder(reg, bird, bin)
  
  ## --- movement events + heading from the regular grid ---
  reg[, `:=`(xp = shift(x), yp = shift(y), bp = shift(bin)), by = bird]
  reg <- reg[(bin - bp) == 1L]
  reg[, disp := sqrt((x-xp)^2 + (y-yp)^2)]
  reg[, heading := atan2(y - yp, x - xp)]
  if (nrow(reg) == 0) return(NULL)
  thr <- median(reg$disp)
  reg[, event := as.integer(disp > thr)]
  
  ## --- matrices ---
  reg[, w := as.integer(factor(bird, levels = birds))]
  m <- length(birds); B <- max(reg$bin) + 1L
  Y <- matrix(0, m, B); Hd <- matrix(NA, m, B)
  PX <- matrix(NA, m, B); PY <- matrix(NA, m, B)
  for (r in seq_len(nrow(reg))) {
    i <- reg$w[r]; b <- reg$bin[r] + 1L
    Y[i,b] <- reg$event[r]; Hd[i,b] <- reg$heading[r]
    PX[i,b] <- reg$x[r];    PY[i,b] <- reg$y[r]
  }
  rf <- sapply(1:m, function(i){bb <- reg[w==i, bin]; if(length(bb)) min(bb)+1L else 1L})
  rl <- sapply(1:m, function(i){bb <- reg[w==i, bin]; if(length(bb)) max(bb)+1L else 0L})
  
  ## --- kernels ---
  L <- CFG$L; beta <- CFG$beta; kappa <- CFG$kappa
  sigma_s <- CFG$sigma_s; d_max <- CFG$d_max
  g_t <- exp(-beta*(1:L)); g_t <- g_t/sum(g_t); I0 <- besselI(kappa, 0)
  
  design <- function(i) {
    thi <- Hd[i,]; pxi <- PX[i,]; pyi <- PY[i,]
    cols <- list(); ci <- 1
    for (j in 1:m) { if (j==i) next
      x <- numeric(B)
      for (l in 1:L) {
        yj <- numeric(B); yj[(l+1):B] <- Y[j,1:(B-l)]
        thj <- rep(NA,B); thj[(l+1):B] <- Hd[j,1:(B-l)]
        pxj <- rep(NA,B); pxj[(l+1):B] <- PX[j,1:(B-l)]
        pyj <- rep(NA,B); pyj[(l+1):B] <- PY[j,1:(B-l)]
        wH <- exp(kappa*cos(thi-thj))/I0
        d2 <- (pxi-pxj)^2 + (pyi-pyj)^2
        wS <- exp(-d2/(2*sigma_s^2)) * (d2 <= d_max^2)
        w <- wH*wS; w[is.na(w)] <- 0
        x <- x + g_t[l]*yj*w
      }
      cols[[ci]] <- x; ci <- ci+1
    }
    do.call(cbind, cols)
  }
  
  alpha <- matrix(0, m, m); se_mat <- matrix(0, m, m); inc_var <- numeric(m)
  for (i in 1:m) {
    if (rl[i] < rf[i]) next
    at <- rf[i]:rl[i]; y_w <- Y[i, at]
    Xi <- design(i)[at, , drop = FALSE]
    if (sum(y_w) < CFG$min_events || all(Xi == 0)) next
    src <- setdiff(1:m, i); nz <- which(colSums(abs(Xi)) > 0)
    if (length(nz) == 0) next
    S_i <- Xi[, nz, drop = FALSE]
    fit <- tryCatch(fit_target_path(S_i, y_w), error = function(e) NULL)
    if (is.null(fit)) next
    for (k in seq_along(nz)) alpha[i, src[nz[k]]] <- fit$alpha[k]
    Vaa <- target_cov(S_i, y_w, fit$mu, fit$alpha)
    se_act <- sqrt(pmax(diag(Vaa), 0)); act <- which(fit$alpha > 1e-6)
    for (kk in seq_along(act)) if (kk <= length(se_act))
      se_mat[i, src[nz[act[kk]]]] <- se_act[kk]
    inc_var[i] <- sum(Vaa)
  }
  
  ## --- directional validation: when i follows j, does i move along j's
  ##     recent heading? (pigeons: expect alignment WELL below baseline) ---
  edges <- which(alpha > 1e-6, arr.ind = TRUE); rows <- list()
  for (e in seq_len(nrow(edges))) {
    i <- edges[e,1]; j <- edges[e,2]
    d2 <- (PX[i,]-PX[j,])^2 + (PY[i,]-PY[j,])^2
    cand <- which(is.finite(Hd[i,]) & Y[i,]==1 & is.finite(d2) & d2 <= d_max^2)
    cand <- cand[cand > 1]; if (!length(cand)) next
    hjp <- Hd[j, cand-1]; ok <- is.finite(hjp); cand <- cand[ok]; hjp <- hjp[ok]
    if (length(cand) < 5) next
    dth <- ang_diff(Hd[i,cand], hjp) * 180/pi
    rows[[length(rows)+1]] <- data.table(target=i, source=j,
                                         alpha=alpha[i,j], n=length(cand),
                                         dtheta=mean(dth))
  }
  val <- if (length(rows)) rbindlist(rows) else NULL
  
  list(birds = birds, alpha = alpha, se_mat = se_mat, inc_var = inc_var,
       out = colSums(alpha), `in` = rowSums(alpha),
       net = colSums(alpha) - rowSums(alpha), val = val,
       rho = max(Mod(eigen(alpha, only.values = TRUE)$values)))
}

## =====================================================================
## 4. RUN ALL FLIGHTS + POOL LEADERSHIP ACROSS RELEASES
##    Nagy's positive-control signature: leadership ranks are CONSISTENT
##    across flights. We pool per-flight net leadership into a per-bird
##    average rank and report cross-flight rank consistency.
## =====================================================================
flights <- split(all_raw, all_raw$frkey)
res <- list(); val_all <- list()
for (fk in names(flights)) {
  r <- tryCatch(analyze_flight(copy(flights[[fk]])), error = function(e){
    message("flight ", fk, ": ", conditionMessage(e)); NULL })
  if (is.null(r)) next
  res[[fk]] <- r
  net <- data.table(frkey = fk, flock = sub("_.*","",fk),
                    bird = r$birds, net = r$net,
                    rank = rank(-r$net))                ## 1 = top leader
  res[[fk]]$net_dt <- net
  if (!is.null(r$val)) { v <- copy(r$val); v[, frkey := fk]; val_all[[fk]] <- v }
  cat(sprintf("flight %-10s | birds %d | rho %.3f | edges %s\n",
              fk, length(r$birds), r$rho,
              if (is.null(r$val)) 0 else nrow(r$val)))
}

net_long <- rbindlist(lapply(res, `[[`, "net_dt"))
## per-bird mean net leadership and mean rank, within each flock
pooled <- net_long[, .(n_flights = .N,
                       mean_net = mean(net),
                       mean_rank = mean(rank)),
                   by = .(flock, bird)][order(flock, mean_rank)]
cat("\n=== pooled leadership by flock (mean_rank: 1 = leader) ===\n")
print(pooled)

## cross-flight rank consistency within each flock (Nagy reported r ~ 0.5-0.6).
## Pairwise Spearman between releases on the shared birds:
consistency <- net_long[, {
  w <- dcast(.SD, bird ~ frkey, value.var = "rank")
  M <- as.matrix(w[, -1, with = FALSE])
  if (ncol(M) >= 2) {
    cc <- cor(M, use = "pairwise.complete.obs", method = "spearman")
    .(mean_pairwise_rho = mean(cc[upper.tri(cc)], na.rm = TRUE),
      n_flights = ncol(M))
  } else .(mean_pairwise_rho = NA_real_, n_flights = ncol(M))
}, by = flock]
cat("\n=== within-flock cross-release rank consistency (Spearman) ===\n")
print(consistency)

## directional validation pooled over all flights + a random-pair baseline
val_dt <- if (length(val_all)) rbindlist(val_all) else NULL
if (!is.null(val_dt)) {
  cat("\n=== directional validation (deg; LOWER = follower tracks leader) ===\n")
  cat("strong edges (top quartile alpha) mean |dtheta|:",
      round(val_dt[alpha >= quantile(alpha,.75), mean(dtheta)], 1), "\n")
  cat("all edges median |dtheta|:", round(median(val_dt$dtheta), 1), "\n")
  cat("NOTE: compute a random-pair baseline per flight to compare against;\n",
      "for pigeons we expect strong-edge dtheta well below the baseline.\n")
}

## =====================================================================
## 5. (optional) save the per-flight alpha matrices for figures
## =====================================================================
saveRDS(res, "pigeon_hawkes_results.rds")
cat("\nDone. Inspect `pooled`, `consistency`, and `val_dt`.\n")

