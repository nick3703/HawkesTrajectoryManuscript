#!/usr/bin/env Rscript
## =====================================================================
##  HawkesTrajAnalysis.R
##
##  Regenerates every numbered table and figure in
##  "Hawkes Trajectories: A Framework for Leadership Inference in
##   Coordinated Movement".
##
##  Each block is named for what it produces in the manuscript.
##
##    Tables1to6      Tables 1-6   simulation studies: recovery with fixed
##                                 kernels, estimated kernels, robustness,
##                                 head-to-head against directional delay
##    Table7          Table 7      bin width and working likelihood
##    Table8          Table 8      sandwich coverage under a binary response
##    Sec6Edges       Section 6    threshold-free edge selection, AUC and
##                                 average precision, bootstrap rule
##    Sec6Driver      Section 6    shared group-level driver simulation
##    Sec6Calib       Section 6    calibrated coupling sweep
##    Sec7Fits        Section 7    intensity decomposition, rho(R) distribution
##    Sec7Boot        Section 7    cluster bootstrap intervals
##    Sec7Resid       Section 7    residual diagnostics
##    Sec7Lag         Section 7    2 s oscillation and its shared component
##    Sec7Origin      Section 7    raw-signal spectrum before binning
##    Table9          Table 9      null battery
##    Sec7Blocks      Section 7    block-restricted timing nulls
##    Sec7Thresh      Section 7    null battery across event thresholds
##    Sec7Events      Section 7    alternative event definitions
##    Sec7RhoACF      Section 7    rho(R) against residual structure
##    Table10         Table 10     group-activity decomposition
##    Sec7FlockNull   Section 7    null battery under the group-activity model
##    Figure1         Figure 1     representative flight
##
##  USAGE
##    Rscript HawkesTrajAnalysis.R            # everything
##    Rscript HawkesTrajAnalysis.R Table9     # one block
##    Rscript HawkesTrajAnalysis.R Table9 Figure1
##
##  SEEDS. SEED_BASE is set once below and each block adds a fixed offset, so
##  any block reproduces identically whether run alone or as part of a full
##  run. Setting the seed only once at the top would not have this property,
##  since the state on entering a block would depend on which blocks preceded
##  it.
##
##  Pigeon data (Watts et al. 2016, Dryad doi:10.5061/dryad.508j2) in ./Data.
##  Dependencies: data.table, sf.
## =====================================================================

suppressPackageStartupMessages({ library(data.table); library(sf) })
SEED_BASE <- 20260806
set.seed(SEED_BASE)

ARGS   <- commandArgs(trailingOnly = TRUE)
BLOCKS <- if (length(ARGS)) ARGS else
  c("Tables1to6","Table7","Table8","Sec6Edges","Sec6Driver","Sec6Calib",
    "Sec7Fits","Sec7Boot","Sec7Resid","Sec7Lag","Sec7Origin",
    "Table9","Sec7Blocks","Sec7Thresh","Sec7Events","Sec7RhoACF",
    "Table10","Sec7FlockNull","Figure1")
run <- function(b) b %in% BLOCKS
DATA_DIR <- "Data"; OUT_DIR <- "."
dir.create(file.path(OUT_DIR,"ckpt"), showWarnings=FALSE, recursive=TRUE)

CFG <- list(
  pattern="\\.csv$", utm_epsg=32630,
  speed_kmh_fly=15, min_fly_sec=30,
  bin_sec=0.4, L=8,
  beta=3.0, kappa=6, sigma_d=15, d_max=60,
  thr_q=0.50, min_events=20,
  eps_mu=1e-8, tol_edge=1e-6, min_colmass=1e-3,
  ## simulation
  M=5, B=600, mu=0.50, coh_dist=30, coh_w=0.5,
  NREP=200, NREP_ETA=200,
  ## nulls
  N_SWAP=20, N_PERM=100, R_BOOT=400,
  QUANTILES=c(0.25,0.50,0.75,0.90),
  ## misc
  n_decile=10, nagy_maxlag=8, nagy_hp_win=15, fs=5
)


## =====================================================================
## PART 1  SHARED: kernels, helpers, metrics
## =====================================================================
gt_of <- function(beta,cfg){ g<-exp(-beta*(1:cfg$L)*cfg$bin_sec); g/sum(g) }
wh_of <- function(dth,k)     exp(k*(cos(dth)-1))
wd_of <- function(d2,s,dm)   exp(-d2/(2*s^2))*(d2<=dm^2)
angd  <- function(a,b)       atan2(sin(a-b),cos(a-b))
lagv  <- function(v,l,B){o<-rep(NA_real_,B); o[(l+1):B]<-v[1:(B-l)]; o}
lagi  <- function(v,l,B){o<-rep(0,B);        o[(l+1):B]<-v[1:(B-l)]; o}
g_t   <- gt_of(CFG$beta,CFG)

mcse   <- function(x) sd(x,na.rm=TRUE)/sqrt(sum(is.finite(x)))
mcse_p <- function(p,n) sqrt(p*(1-p)/n)
emp_p  <- function(r,n) if(!length(n)) NA_real_ else (1+sum(n>=r))/(1+length(n))
fisher <- function(p){p<-p[is.finite(p)]
if(!length(p)) NA_real_ else pchisq(-2*sum(log(pmax(p,1e-6))),2*length(p),lower.tail=FALSE)}
med    <- function(x) if(length(x)) median(x,na.rm=TRUE) else NA_real_

net_lead <- function(A){off<-A; diag(off)<-0; colSums(off)-rowSums(off)}
rank_ceiling <- function(nl){o<-order(nl,decreasing=TRUE)
id<-numeric(length(nl)); id[o]<-length(nl):1
suppressWarnings(cor(nl,id,method="spearman"))}
metrics <- function(Ah,At,tol){
  oh<-Ah; ot<-At; diag(oh)<-0; diag(ot)<-0
  pos<-ot[ot>0]; thr<-if(length(pos)) 0.5*min(pos) else tol
  det<-oh>thr; tru<-ot>0
  TP<-sum(det&tru); FP<-sum(det&!tru); FN<-sum(!det&tru)
  nt<-net_lead(At); nh<-net_lead(Ah); top<-which(nt>=max(nt)-1e-12)
  list(sens=if(TP+FN>0)TP/(TP+FN) else NA_real_,
       prec=if(TP+FP>0)TP/(TP+FP) else NA_real_,
       fdr =if(TP+FP>0)FP/(TP+FP) else NA_real_,
       edge_rmse=sqrt(mean((oh-ot)^2)), edge_medae=median(abs(oh-ot)),
       max_alpha=max(oh), null_edge_mean=if(any(!tru)) mean(oh[!tru]) else NA_real_,
       self_bias=mean(diag(Ah)-diag(At)),
       rank_cor=suppressWarnings(cor(nt,nh,method="spearman")),
       top_correct=which.max(nh)%in%top)
}

## core estimator: identity-link Poisson, mu>=eps, alpha>=0, alpha_ii free
fit_row <- function(S,y,cfg,want_se=FALSE,want_fit=FALSE){
  keep<-which(colMeans(S)>cfg$min_colmass); if(!length(keep)) return(NULL)
  A1<-cbind(1,S[,keep,drop=FALSE]); p<-length(keep)
  nll<-function(th){lam<-pmax(as.numeric(A1%*%th),1e-12);-sum(y*log(lam)-lam)}
  gr <-function(th){lam<-pmax(as.numeric(A1%*%th),1e-12)
  -as.numeric(crossprod(A1,y/lam-1))}
  o<-tryCatch(optim(c(max(mean(y),1e-6),rep(.1,p)),nll,gr,method="L-BFGS-B",
                    lower=c(cfg$eps_mu,rep(0,p)),
                    control=list(maxit=400,factr=1e10)),error=function(e)NULL)
  if(is.null(o)) return(NULL)
  a<-numeric(ncol(S)); a[keep]<-o$par[-1]
  out<-list(mu=o$par[1],alpha=a,keep=keep)
  lam<-pmax(as.numeric(A1%*%o$par),1e-12)
  out$loglik<-sum(y*log(lam)-lam); out$npar<-p+1
  if(want_fit) out$lambda<-lam
  if(want_se){ I<-crossprod(A1,A1/lam)
  V<-tryCatch(solve(I),error=function(e)NULL)
  if(!is.null(V)){s<-sqrt(pmax(diag(V),0)); se<-numeric(ncol(S)); se[keep]<-s[-1]
  out$se<-se; out$cond<-tryCatch(kappa(I),error=function(e)NA_real_)}}
  out
}

fit_row_ql <- function(S, y, cfg, hac_lag = cfg$L) {
  keep <- which(colMeans(S) > cfg$min_colmass); if (!length(keep)) return(NULL)
  A1 <- cbind(1, S[, keep, drop = FALSE]); p <- length(keep)
  nll <- function(th){lam <- pmax(as.numeric(A1 %*% th), 1e-12); -sum(y*log(lam)-lam)}
  gr  <- function(th){lam <- pmax(as.numeric(A1 %*% th), 1e-12)
  -as.numeric(crossprod(A1, y/lam - 1))}
  o <- tryCatch(optim(c(max(mean(y),1e-6), rep(.1,p)), nll, gr, method = "L-BFGS-B",
                      lower = c(cfg$eps_mu, rep(0,p)),
                      control = list(maxit = 400, factr = 1e10)), error = function(e) NULL)
  if (is.null(o)) return(NULL)
  th  <- o$par
  lam <- pmax(as.numeric(A1 %*% th), 1e-12)
  r   <- (y - lam) / lam                      # score multiplier
  U   <- A1 * r                               # n x (p+1) score contributions
  
  J <- crossprod(A1, A1 / lam)                # bread
  K_rob <- crossprod(U)                       # meat, independent across bins
  ## HAC meat: Gamma_0 + sum_h w_h (Gamma_h + Gamma_h'), Bartlett weights
  K_hac <- K_rob
  if (hac_lag > 0) for (h in seq_len(hac_lag)) {
    n <- nrow(U); if (h >= n) break
    G <- crossprod(U[(h+1):n, , drop = FALSE], U[1:(n-h), , drop = FALSE])
    K_hac <- K_hac + (1 - h/(hac_lag + 1)) * (G + t(G))
  }
  Ji <- tryCatch(solve(J), error = function(e) NULL)
  se <- function(V) if (is.null(V)) rep(NA_real_, p+1) else sqrt(pmax(diag(V), 0))
  V_naive <- Ji
  V_rob   <- if (is.null(Ji)) NULL else Ji %*% K_rob %*% Ji
  V_hac   <- if (is.null(Ji)) NULL else Ji %*% K_hac %*% Ji
  expand <- function(v){ out <- rep(NA_real_, ncol(S)); out[keep] <- v[-1]; out }
  a <- numeric(ncol(S)); a[keep] <- th[-1]
  list(mu = th[1], alpha = a, keep = keep, lambda = lam,
       se_naive = expand(se(V_naive)), se_rob = expand(se(V_rob)),
       se_hac = expand(se(V_hac)),
       loglik = sum(y*log(lam) - lam), npar = p + 1,
       dispersion = sum((y - lam)^2 / lam) / max(length(y) - p - 1, 1))
}
SD <- list(
  NREP     = 200,     # replicates per arm
  N_PERM   = 30,      # timing-null replicates per bird (floor p = 1/31 = 0.032)
  H        = 10,      # window for the group-activity term, bins
  drive_sd = 0.6,     # log-scale amplitude of the common driver
  drive_ac = 0.85,    # AR(1) coefficient of the driver; 0.85 at 0.4 s bins
  # gives a correlation time of about 2.5 s, comparable to
  # the oscillation seen in the pigeon residuals
  alpha_x  = 0.30     # off-diagonal excitation in the arms that have coupling
)

## common driver: a positive, autocorrelated multiplier shared by all birds.
## Predictable by construction, since bird i's intensity in bin b uses the
## driver value at b, which does not depend on any bird's events.
make_driver <- function(B, cfg) {
  e <- numeric(B); e[1] <- rnorm(1, 0, SD$drive_sd)
  for (b in 2:B) e[b] <- SD$drive_ac*e[b-1] + rnorm(1, 0, SD$drive_sd*sqrt(1-SD$drive_ac^2))
  exp(e - max(e))      # scaled to at most one so the baseline rate is not inflated
}

## generate events. `drive` multiplies the background; `A` supplies excitation.
simulate_arm <- function(W, A, drv, cfg) {
  M <- cfg$M; B <- cfg$B; Y <- matrix(0, M, B)
  for (b in (cfg$L+1):B) for (i in seq_len(M)) {
    lam <- cfg$mu*cfg$bin_sec*drv[b]        # driver enters through the background
    for (j in seq_len(M)) for (l in seq_len(cfg$L)) {
      wv <- W[[i]][[j]][[l]][b]
      if (is.finite(wv) && wv != 0) lam <- lam + A[i,j]*wv*Y[j,b-l]
    }
    Y[i,b] <- rbinom(1, 1, min(max(lam,0), 1))
  }
  Y
}

## group-activity regressor, as in the application
z_of <- function(Y, i, H, B) {
  m <- nrow(Y); others <- setdiff(seq_len(m), i)
  acc <- numeric(B)
  for (j in others) for (h in seq_len(H)) acc <- acc + lagi(Y[j,], h, B)
  acc/(length(others)*H)
}

## fit one replicate, optionally with the group-activity term
fit_arm <- function(W, Y, cfg, use_gamma, Yover = NULL) {
  Yf <- if (is.null(Yover)) Y else Yover
  m <- cfg$M; A <- matrix(0, m, m); at <- (cfg$L+1):cfg$B
  for (i in seq_len(m)) {
    S <- design_sim(W, Yf, i, cfg)
    X <- if (use_gamma) cbind(S, z_of(Yf, i, SD$H, cfg$B)) else S
    f <- fit_row(X[at,,drop=FALSE], Yf[i,at], cfg)
    if (is.null(f)) return(NULL)
    A[i,] <- f$alpha[seq_len(m)]
  }
  off <- A; diag(off) <- 0
  list(off=off, out=colSums(off), inc=rowSums(off))
}

## timing nulls on source b, as in the application
perm_row <- function(v) { v[sample(length(v))] }
circ_row <- function(v, L) {
  n <- length(v); if (n < 2*L+3) return(NULL)
  k <- sample(seq(L+1, n-L-1), 1); v[((seq_len(n)-1+k) %% n) + 1]
}

## =====================================================================
## PART 2  SIMULATION MACHINERY
## =====================================================================
make_A <- function(sc,M,self=0.40){
  A<-matrix(0,M,M); diag(A)<-self
  switch(sc,
         "single"  ={A[2:M,1]<-0.30},
         "chain"   ={for(i in 2:M) A[i,i-1]<-0.35},
         "coleader"={A[3:M,1]<-0.25; A[3:M,2]<-0.25},
         "subgroup"={A[2,1]<-A[3,1]<-0.30; A[5,4]<-0.30},
         stop("bad scenario"))
  stopifnot(max(rowSums(A))<1); A
}
SCEN <- c("single","chain","coleader","subgroup")
stopifnot(make_A("chain",CFG$M)[2,1]>0, make_A("chain",CFG$M)[1,2]==0)  # rows=targets

simulate_marks <- function(cfg,align=0.6){
  M<-cfg$M;B<-cfg$B
  Hd<-PX<-PY<-matrix(NA_real_,M,B)
  common<-cumsum(rnorm(B,0,0.05))
  hd<-rnorm(M,0,0.3); x<-rnorm(M,0,10); y<-rnorm(M,0,10); spd<-15*cfg$bin_sec
  for(b in seq_len(B)){
    hd<-align*(common[b]+rnorm(M,0,0.15))+(1-align)*(hd+rnorm(M,0,0.25))
    cx<-mean(x);cy<-mean(y); dd<-sqrt((cx-x)^2+(cy-y)^2); br<-atan2(cy-y,cx-x)
    p<-pmin(dd/cfg$coh_dist,1)*cfg$coh_w
    hd<-atan2((1-p)*sin(hd)+p*sin(br),(1-p)*cos(hd)+p*cos(br))
    x<-x+spd*cos(hd); y<-y+spd*sin(hd)
    Hd[,b]<-hd; PX[,b]<-x; PY[,b]<-y
  }
  list(Hd=Hd,PX=PX,PY=PY)
}
build_W_sim <- function(mk,eta,cfg){
  M<-cfg$M;B<-cfg$B; g<-gt_of(eta$beta,cfg)
  tg<-lapply(seq_len(M),function(i)
    list(th=lagv(mk$Hd[i,],1,B),px=lagv(mk$PX[i,],1,B),py=lagv(mk$PY[i,],1,B)))
  W<-vector("list",M)
  for(i in seq_len(M)){W[[i]]<-vector("list",M)
  for(j in seq_len(M)){wl<-vector("list",cfg$L)
  for(l in seq_len(cfg$L)){
    hj<-lagv(mk$Hd[j,],l,B);pxj<-lagv(mk$PX[j,],l,B);pyj<-lagv(mk$PY[j,],l,B)
    d2<-(tg[[i]]$px-pxj)^2+(tg[[i]]$py-pyj)^2
    w<-wh_of(angd(tg[[i]]$th,hj),eta$kappa)*wd_of(d2,eta$sigma_d,cfg$d_max)
    w[!is.finite(w)]<-0; wl[[l]]<-g[l]*w}
  W[[i]][[j]]<-wl}}
  W
}
simulate_Y <- function(W,A,mu,cfg,family="poisson"){
  M<-cfg$M;B<-cfg$B; Y<-matrix(0,M,B)
  for(b in (cfg$L+1):B) for(i in seq_len(M)){
    lam<-mu*cfg$bin_sec
    for(j in seq_len(M)) for(l in seq_len(cfg$L)){
      wv<-W[[i]][[j]][[l]][b]
      if(is.finite(wv)&&wv!=0) lam<-lam+A[i,j]*wv*Y[j,b-l]}
    Y[i,b]<-if(family=="poisson") rpois(1,max(lam,0)) else rbinom(1,1,min(max(lam,0),1))
  }
  Y
}
design_sim <- function(W,Y,i,cfg){
  S<-matrix(0,cfg$B,cfg$M)
  for(j in seq_len(cfg$M)) for(l in seq_len(cfg$L))
    S[,j]<-S[,j]+W[[i]][[j]][[l]]*lagi(Y[j,],l,cfg$B)
  S
}
fit_sim <- function(W,Y,cfg,want_se=FALSE){
  A<-matrix(0,cfg$M,cfg$M); SE<-matrix(NA_real_,cfg$M,cfg$M)
  at<-(cfg$L+1):cfg$B
  for(i in seq_len(cfg$M)){
    f<-fit_row(design_sim(W,Y,i,cfg)[at,,drop=FALSE],Y[i,at],cfg,want_se)
    if(is.null(f)) return(NULL)
    A[i,]<-f$alpha; if(want_se&&!is.null(f$se)) SE[i,]<-f$se}
  list(A=A,SE=SE)
}
delay_leader <- function(mk,cfg){
  M<-nrow(mk$Hd);B<-ncol(mk$Hd)
  UX<-cos(mk$Hd);UY<-sin(mk$Hd)
  hp<-function(v){n<-length(v);k<-floor(cfg$nagy_hp_win/2);ma<-rep(NA_real_,n)
  cs<-c(0,cumsum(ifelse(is.na(v),0,v)));ct<-c(0,cumsum(!is.na(v)))
  for(t in seq_len(n)){lo<-max(1,t-k);hi<-min(n,t+k);cn<-ct[hi+1]-ct[lo]
  ma[t]<-if(cn>0)(cs[hi+1]-cs[lo])/cn else NA_real_}; v-ma}
  FX<-t(apply(UX,1,hp));FY<-t(apply(UY,1,hp))
  nr<-sqrt(FX^2+FY^2); nr[!is.finite(nr)|nr<1e-12]<-NA; FX<-FX/nr;FY<-FY/nr
  sc<-numeric(M)
  for(i in seq_len(M)) for(j in seq_len(M)){ if(i==j) next
    best<--Inf;bt<-0
    for(tau in (-cfg$nagy_maxlag):cfg$nagy_maxlag){
      if(tau>=0){a<-1:(B-tau);b2<-(1+tau):B} else {s<--tau;a<-(1+s):B;b2<-1:(B-s)}
      cc<-mean(FX[i,a]*FX[j,b2]+FY[i,a]*FY[j,b2],na.rm=TRUE)
      if(is.finite(cc)&&cc>best){best<-cc;bt<-tau}}
    if(bt>0) sc[i]<-sc[i]+1 else if(bt<0) sc[j]<-sc[j]+1}
  sc
}
simulate_follow <- function(cfg,fid,sdn=0.35){
  M<-cfg$M;B<-cfg$B; Hd<-PX<-PY<-matrix(NA_real_,M,B)
  hd<-rnorm(M,0,.3);x<-rnorm(M,0,10);y<-rnorm(M,0,10)
  base<-15*cfg$bin_sec; spd<-rep(base,M); prev<-base
  for(b in seq_len(B)){
    lead<-hd[1]+rnorm(1,0,0.30); spd[1]<-base*exp(rnorm(1,0,sdn))
    for(i in 2:M){
      hd[i]<-if(runif(1)<fid) lead+rnorm(1,0,0.15) else hd[i]+rnorm(1,0,0.30)
      spd[i]<-exp(fid*log(prev)+(1-fid)*log(base)+rnorm(1,0,sdn))}
    hd[1]<-lead
    cx<-mean(x);cy<-mean(y); dd<-sqrt((cx-x)^2+(cy-y)^2); br<-atan2(cy-y,cx-x)
    pp<-pmin(dd/cfg$coh_dist,1)*cfg$coh_w
    hd<-atan2((1-pp)*sin(hd)+pp*sin(br),(1-pp)*cos(hd)+pp*cos(br))
    x<-x+spd*cos(hd); y<-y+spd*sin(hd)
    Hd[,b]<-hd;PX[,b]<-x;PY[,b]<-y; prev<-spd[1]}
  disp<-cbind(0,t(sqrt(diff(t(PX))^2+diff(t(PY))^2)))
  list(mk=list(Hd=Hd,PX=PX,PY=PY),
       Y=(disp>matrix(apply(disp,1,median),M,B))*1)
}

## =====================================================================
## PART 3  PIGEON DATA PIPELINE
## =====================================================================
load_pigeons <- function(dir,cfg){
  files<-list.files(dir,pattern=cfg$pattern,full.names=TRUE,recursive=TRUE)
  if(!length(files)) stop("No CSVs under ",dir)
  one<-function(p){
    part<-strsplit(sub("\\.csv$","",basename(p)),"_")[[1]]
    d<-fread(p); setnames(d,trimws(names(d))); d<-d[VALID=="FIXED"]
    if(!nrow(d)) return(NULL)
    d[,time:=as.POSIXct(paste(`UTC DATE`,`UTC TIME`),tz="UTC",
                        format="%Y/%m/%d %H:%M:%S")+MS/1000]
    d[,lat:=fifelse(`N/S`=="S",-1,1)*LATITUDE]
    d[,lon:=fifelse(`E/W`=="W",-1,1)*LONGITUDE]
    d[,.(flock=part[2],bird=part[3],release=part[4],time,lon,lat,speed=SPEED)]}
  A<-rbindlist(lapply(files,one),use.names=TRUE)
  setorder(A,flock,release,bird,time); A[,frkey:=paste(flock,release,sep="_")]
  xy<-sf::sf_project("EPSG:4326",paste0("EPSG:",cfg$utm_epsg),as.matrix(A[,.(lon,lat)]))
  A[,`:=`(x=xy[,1],y=xy[,2])]; A[]
}

## grid one flight; events defined at quantile q of `evtype`
prep_flight <- function(D,cfg,q=cfg$thr_q,evtype="displacement"){
  birds<-sort(unique(D$bird)); if(length(birds)<2) return(NULL)
  D[,fly:=speed>=cfg$speed_kmh_fly]; flo<-D[fly==TRUE]; if(!nrow(flo)) return(NULL)
  tlo<-min(flo$time);thi<-max(flo$time)
  if(as.numeric(thi-tlo,units="secs")<cfg$min_fly_sec) return(NULL)
  D<-D[time>=tlo&time<=thi]; t0<-as.numeric(min(D$time))
  D[,bin:=as.integer(floor((as.numeric(time)-t0)/cfg$bin_sec))]
  reg<-D[,.(x=mean(x),y=mean(y)),by=.(bird,bin)]; setorder(reg,bird,bin)
  reg[,`:=`(xp=shift(x),yp=shift(y),bp=shift(bin)),by=bird]
  reg<-reg[(bin-bp)==1L]; if(!nrow(reg)) return(NULL)
  reg[,disp:=sqrt((x-xp)^2+(y-yp)^2)]
  reg[,heading:=atan2(y-yp,x-xp)]
  reg[,dspeed:=abs(disp-shift(disp)),by=bird]
  reg[,dhead :=abs(angd(heading,shift(heading))),by=bird]
  v<-switch(evtype,displacement="disp",acceleration="dspeed",turning="dhead")
  reg[,.val:=get(v)]; reg<-reg[is.finite(.val)]; if(!nrow(reg)) return(NULL)
  reg[,event:=as.integer(.val>quantile(.val,q,names=FALSE)),by=bird]
  m<-length(birds); B<-max(reg$bin)+1L
  Y<-OBS<-matrix(0L,m,B); Hd<-PX<-PY<-matrix(NA_real_,m,B)
  wv<-as.integer(factor(reg$bird,levels=birds))
  for(r in seq_len(nrow(reg))){i<-wv[r];b<-reg$bin[r]+1L; if(b<1||b>B) next
  Y[i,b]<-reg$event[r];OBS[i,b]<-1L
  Hd[i,b]<-reg$heading[r];PX[i,b]<-reg$x[r];PY[i,b]<-reg$y[r]}
  list(birds=birds,m=m,B=B,Y=Y,OBS=OBS,Hd=Hd,PX=PX,PY=PY,reg=reg,
       centroid=reg[,.(cx=mean(x),cy=mean(y)),by=bin])
}
mats_from_reg <- function(reg,birds,B){
  m<-length(birds); Y<-OBS<-matrix(0L,m,B); Hd<-PX<-PY<-matrix(NA_real_,m,B)
  wv<-as.integer(factor(reg$bird,levels=birds))
  for(r in seq_len(nrow(reg))){i<-wv[r];b<-reg$bin[r]+1L; if(b<1||b>B) next
  Y[i,b]<-reg$event[r];OBS[i,b]<-1L
  Hd[i,b]<-reg$heading[r];PX[i,b]<-reg$x[r];PY[i,b]<-reg$y[r]}
  list(m=m,B=B,Y=Y,OBS=OBS,Hd=Hd,PX=PX,PY=PY)
}
build_W <- function(M,cfg,predictable=TRUE){
  m<-M$m;B<-M$B
  tg<-lapply(seq_len(m),function(i) if(predictable)
    list(th=lagv(M$Hd[i,],1,B),px=lagv(M$PX[i,],1,B),py=lagv(M$PY[i,],1,B))
    else list(th=M$Hd[i,],px=M$PX[i,],py=M$PY[i,]))
  W<-vector("list",m); rows<-vector("list",m)
  for(i in seq_len(m)){W[[i]]<-vector("list",m)
  for(j in seq_len(m)){wl<-vector("list",cfg$L)
  for(l in seq_len(cfg$L)){
    hj<-lagv(M$Hd[j,],l,B);pxj<-lagv(M$PX[j,],l,B);pyj<-lagv(M$PY[j,],l,B)
    d2<-(tg[[i]]$px-pxj)^2+(tg[[i]]$py-pyj)^2
    w<-wh_of(angd(tg[[i]]$th,hj),cfg$kappa)*wd_of(d2,cfg$sigma_d,cfg$d_max)
    w[!is.finite(w)]<-0; wl[[l]]<-g_t[l]*w}
  W[[i]][[j]]<-wl}
  at<-which(M$OBS[i,]==1L); rows[[i]]<-at[at>cfg$L]}
  list(W=W,rows=rows,m=m,B=B)
}
col_of <- function(DZ,i,j,Yj,cfg){x<-numeric(DZ$B)
for(l in seq_len(cfg$L)) x<-x+DZ$W[[i]][[j]][[l]]*lagi(Yj,l,DZ$B); x}
design_i <- function(DZ,Y,i,cfg){
  S<-matrix(0,DZ$B,DZ$m); for(j in seq_len(DZ$m)) S[,j]<-col_of(DZ,i,j,Y[j,],cfg); S}
S0_i <- function(DZ,Y,i,cfg){
  S<-matrix(0,DZ$B,DZ$m)
  for(j in seq_len(DZ$m)) for(l in seq_len(cfg$L)) S[,j]<-S[,j]+g_t[l]*lagi(Y[j,],l,DZ$B); S}

fit_flight <- function(DZ,M,cfg,targets=seq_len(DZ$m),Sover=NULL,
                       want_se=FALSE,want_fit=FALSE,want_Ew=FALSE){
  m<-DZ$m; alpha<-matrix(0,m,m); SE<-Ew<-matrix(NA_real_,m,m)
  mu<-rep(NA_real_,m); cond<-rep(NA_real_,m); ll<-0; np<-0
  sh<-matrix(NA_real_,m,3); FIT<-list()
  for(i in targets){
    at<-DZ$rows[[i]]; if(length(at)<2) next
    y<-M$Y[i,at]; if(sum(y)<cfg$min_events) next
    S<-if(is.null(Sover)) design_i(DZ,M$Y,i,cfg) else Sover[[i]]
    Si<-S[at,,drop=FALSE]
    f<-fit_row(Si,y,cfg,want_se,want_fit); if(is.null(f)) next
    alpha[i,]<-f$alpha; mu[i]<-f$mu; ll<-ll+f$loglik; np<-np+f$npar
    if(want_se&&!is.null(f$se)){SE[i,]<-f$se; cond[i]<-f$cond}
    if(want_fit) FIT[[length(FIT)+1L]]<-data.table(bird=i,bin=at,
                                                   t_rel=(at-min(at))/max(1,(max(at)-min(at))),y=y,lambda=f$lambda)
    if(want_Ew){S0<-S0_i(DZ,M$Y,i,cfg)[at,,drop=FALSE]
    for(j in seq_len(m)){m0<-mean(S0[,j]); Ew[i,j]<-if(m0>1e-12) mean(Si[,j])/m0 else NA_real_}}
    con<-sweep(Si,2,f$alpha,`*`); tot<-mean(f$mu+rowSums(con))
    if(is.finite(tot)&&tot>0) sh[i,]<-c(f$mu/tot,mean(con[,i])/tot,
                                        (mean(rowSums(con))-mean(con[,i]))/tot)
  }
  off<-alpha; diag(off)<-0
  R<-alpha*ifelse(is.na(Ew),0,Ew)
  list(alpha=alpha,off=off,mu=mu,SE=SE,cond=cond,Ew=Ew,
       out=colSums(off),inc=rowSums(off),net=colSums(off)-rowSums(off),
       self=diag(alpha),shares=colMeans(sh,na.rm=TRUE),loglik=ll,npar=np,
       n_edges=sum(off>cfg$tol_edge),
       rho_A=max(Mod(eigen(alpha,only.values=TRUE)$values)),
       rho_R=max(Mod(eigen(R,only.values=TRUE)$values)),
       rho_off=max(Mod(eigen(off,only.values=TRUE)$values)),
       fitted=if(want_fit&&length(FIT)) rbindlist(FIT) else NULL)
}

## null generators
perm_ev  <- function(Yr,OBSr){o<-which(OBSr==1L); if(length(o)<2) return(Yr)
Yr[o]<-sample(Yr[o]); Yr}
shift_ev <- function(Yr,OBSr,L){o<-which(OBSr==1L);n<-length(o)
if(n<2*L+3) return(NULL); k<-sample(seq(L+1,n-L-1),1)
Yr[o]<-Yr[o][((seq_len(n)-1+k)%%n)+1]; Yr}
shift_tr <- function(v,OBSr,k){o<-which(OBSr==1L); if(length(o)<3) return(v)
v[o]<-v[o][((seq_along(o)-1+k)%%length(o))+1]; v}
donor_reg <- function(G,donor_full,bird_id,cfg){
  dr<-copy(donor_full[bird==bird_id]); if(!nrow(dr)) return(NULL)
  dr[,bin:=bin-min(bin)]; dr<-dr[bin<G$B]
  if(nrow(dr)<cfg$min_events*2) return(NULL)
  dc<-dr[,.(mx=mean(x),my=mean(y))]
  dr<-merge(dr,G$centroid,by="bin",all.x=TRUE)
  dr[,x:=x-dc$mx+fifelse(is.finite(cx),cx,dc$mx)]
  dr[,y:=y-dc$my+fifelse(is.finite(cy),cy,dc$my)]
  dr[,`:=`(cx=NULL,cy=NULL)]; dr[,bird:=bird_id]
  dr[,.(bird,bin,x,y,disp,heading,event)]}
out_null <- function(DZ,M,b,type,cfg){
  if(type=="permmark"){
    Mm<-M; k<-sample(seq_len(max(sum(M$OBS[b,]),2)),1)
    Mm$Y[b,]<-perm_ev(M$Y[b,],M$OBS[b,])
    Mm$Hd[b,]<-shift_tr(M$Hd[b,],M$OBS[b,],k)
    Mm$PX[b,]<-shift_tr(M$PX[b,],M$OBS[b,],k)
    Mm$PY[b,]<-shift_tr(M$PY[b,],M$OBS[b,],k)
    return(fit_flight(build_W(Mm,cfg),Mm,cfg,targets=setdiff(seq_len(DZ$m),b))$out[b])}
  Yb<-if(type=="perm") perm_ev(M$Y[b,],M$OBS[b,]) else shift_ev(M$Y[b,],M$OBS[b,],cfg$L)
  if(is.null(Yb)) return(NA_real_)
  tg<-setdiff(seq_len(DZ$m),b)
  So<-lapply(seq_len(DZ$m),function(i) design_i(DZ,M$Y,i,cfg))
  for(i in tg) So[[i]][,b]<-col_of(DZ,i,b,Yb,cfg)
  fit_flight(DZ,M,cfg,targets=tg,Sover=So)$out[b]}

PIG <- NULL
need_pigeons <- function(){
  if(is.null(PIG)) PIG <<- load_pigeons(DATA_DIR,CFG)
  invisible(PIG)}
prep_all <- function(q=CFG$thr_q,evtype="displacement"){
  need_pigeons()
  P<-lapply(unique(PIG$frkey),function(fk)
    tryCatch(prep_flight(PIG[frkey==fk],CFG,q,evtype),error=function(e)NULL))
  names(P)<-unique(PIG$frkey); Filter(Negate(is.null),P)}
flock_of <- function(fk) sub("_.*","",fk)

## =====================================================================
## PART 4  ANALYSES
## =====================================================================

## ---------------------------------------------------------------------
#  Generates Tables 1-8 (Simulation Studies 1-5)

## --------------------------------------------------------------------
##  Tables 1-6: simulation studies, recovery and robustness
## --------------------------------------------------------------------
if (run("Tables1to6")) {
  message("\n### SIMULATION STUDY ###")
  set.seed(SEED_BASE + 1000)
  eta0<-list(beta=CFG$beta,sigma_d=CFG$sigma_d,kappa=CFG$kappa)
  
  # Generates Table 1 (true and recovered A) and Table 2 (recovery metrics)
  S1<-list()
  for(sc in SCEN){
    At<-make_A(sc,CFG$M); keep<-0; dfit<-0; dev<-0; Ms<-list(); Ah<-list(); COV<-list()
    for(r in seq_len(CFG$NREP)){
      mk<-simulate_marks(CFG); W<-build_W_sim(mk,eta0,CFG)
      Y<-simulate_Y(W,At,CFG$mu,CFG,"poisson")
      if(min(rowSums(Y))<10){dev<-dev+1;next}
      f<-fit_sim(W,Y,CFG,want_se=TRUE); if(is.null(f)){dfit<-dfit+1;next}
      keep<-keep+1; Ms[[keep]]<-metrics(f$A,At,CFG$tol_edge); Ah[[keep]]<-f$A
      idx<-which(At>0,arr.ind=TRUE)
      COV[[keep]]<-vapply(seq_len(nrow(idx)),function(k){
        i<-idx[k,1];j<-idx[k,2];s<-f$SE[i,j]
        if(!is.finite(s)||s<=0) return(NA)
        abs(f$A[i,j]-At[i,j])<=1.96*s},logical(1))}
    if(!keep) next
    g<-function(f) vapply(Ms,function(z) as.numeric(z[[f]]),numeric(1))
    tc<-mean(g("top_correct")); cv<-unlist(COV)
    cat("\n",sc,"true A:\n"); print(round(At,3))
    cat("mean estimated A over",keep,"replicates:\n"); print(round(Reduce(`+`,Ah)/keep,3))
    S1[[sc]]<-data.table(scenario=sc,replicates=keep,discarded_fit=dfit,
                         discarded_lowevent=dev,
                         sensitivity=round(mean(g("sens"),na.rm=TRUE),3),sens_mcse=round(mcse(g("sens")),4),
                         precision=round(mean(g("prec"),na.rm=TRUE),3),FDR=round(mean(g("fdr"),na.rm=TRUE),3),
                         edge_RMSE=round(mean(g("edge_rmse")),4),edge_RMSE_mcse=round(mcse(g("edge_rmse")),5),
                         null_edge_mean=round(mean(g("null_edge_mean"),na.rm=TRUE),4),
                         self_bias=round(mean(g("self_bias")),4),
                         rank_cor=round(mean(g("rank_cor"),na.rm=TRUE),3),
                         rank_max=round(rank_ceiling(net_lead(At)),3),
                         top_leader=round(tc,3),top_mcse=round(mcse_p(tc,keep),4),
                         coverage=round(mean(cv,na.rm=TRUE),3),
                         cov_mcse=round(mcse_p(mean(cv,na.rm=TRUE),sum(is.finite(cv))),4))}
  T12<-rbindlist(S1); print(T12); fwrite(T12,file.path(OUT_DIR,"table_sim_recovery.csv"))
  
  # Generates Tables 3-4 (Study 2, estimated kernels: parameters and recovery of A)
  prof_fit <- function(mk,Y,cfg){
    prof<-function(eta){W<-build_W_sim(mk,eta,cfg)
    f<-tryCatch(fit_sim(W,Y,cfg),error=function(e)NULL)
    if(is.null(f)) return(-Inf)
    at<-(cfg$L+1):cfg$B; tot<-0
    for(i in seq_len(cfg$M)){
      r<-fit_row(design_sim(W,Y,i,cfg)[at,,drop=FALSE],Y[i,at],cfg)
      if(is.null(r)) return(-Inf); tot<-tot+r$loglik}
    tot}
    best<-list(ll=-Inf,eta=NULL)
    for(b in c(1.5,3,6,10)) for(s in c(4,8,15,30)) for(k in c(2,5,10,20)){
      e<-list(beta=b,sigma_d=s,kappa=k); ll<-prof(e)
      if(is.finite(ll)&&ll>best$ll) best<-list(ll=ll,eta=e)}
    if(is.null(best$eta)) return(NULL)
    lo<-log(c(0.2,2,0.2)); hi<-log(c(30,cfg$d_max,60))
    po<-tryCatch(optim(log(unlist(best$eta[c("beta","sigma_d","kappa")])),
                       function(lp){if(any(lp<lo)||any(lp>hi)) return(1e12)
                         -prof(list(beta=exp(lp[1]),sigma_d=exp(lp[2]),kappa=exp(lp[3])))},
                       method="Nelder-Mead",control=list(maxit=80,reltol=1e-6)),error=function(e)NULL)
    if(!is.null(po)&&is.finite(po$value)&&-po$value>best$ll)
      best<-list(ll=-po$value,eta=list(beta=exp(po$par[1]),sigma_d=exp(po$par[2]),
                                       kappa=exp(po$par[3])))
    f<-fit_sim(build_W_sim(mk,best$eta,cfg),Y,cfg); if(is.null(f)) return(NULL)
    c(f,list(eta=best$eta))}
  S2<-list()
  for(rg in c("tight","dispersed")){
    cfgr<-modifyList(CFG,if(rg=="tight") list() else list(coh_dist=4*CFG$coh_dist,coh_w=0.12))
    for(sc in c("single","chain")){
      At<-make_A(sc,cfgr$M); keep<-0;disc<-0
      E<-matrix(NA_real_,cfgr$NREP_ETA,3); Ms<-list(); dsep<-c()
      for(r in seq_len(cfgr$NREP_ETA)){
        mk<-simulate_marks(cfgr); W<-build_W_sim(mk,eta0,cfgr)
        Y<-simulate_Y(W,At,cfgr$mu,cfgr,"poisson")
        if(min(rowSums(Y))<10){disc<-disc+1;next}
        dsep<-c(dsep,median(vapply(seq(cfgr$L+1,cfgr$B,by=25),function(b)
          median(dist(cbind(mk$PX[,b],mk$PY[,b]))),numeric(1))))
        f<-prof_fit(mk,Y,cfgr); if(is.null(f)){disc<-disc+1;next}
        keep<-keep+1; E[keep,]<-unlist(f$eta[c("beta","sigma_d","kappa")])
        Ms[[keep]]<-metrics(f$A,At,cfgr$tol_edge)}
      if(!keep) next
      E<-E[seq_len(keep),,drop=FALSE]
      g<-function(f) vapply(Ms,function(z) as.numeric(z[[f]]),numeric(1))
      tru<-c(cfgr$beta,cfgr$sigma_d,cfgr$kappa)
      ab<-c(mean(E[,1]<=0.21|E[,1]>=29.9),mean(E[,2]<=2.05|E[,2]>=cfgr$d_max*.995),
            mean(E[,3]<=0.21|E[,3]>=59.5))
      for(k in 1:3) S2[[paste(rg,sc,k)]]<-data.table(regime=rg,sep_m=round(mean(dsep),1),
                                                     scenario=sc,parameter=c("beta","sigma_d","kappa")[k],truth=tru[k],
                                                     median_est=round(median(E[,k]),2),
                                                     IQR_lo=round(quantile(E[,k],.25,names=FALSE),2),
                                                     IQR_hi=round(quantile(E[,k],.75,names=FALSE),2),
                                                     mean_est=round(mean(E[,k]),2),med_bias=round(median(E[,k])-tru[k],3),
                                                     MCSE=round(mcse(E[,k]),3),pct_at_bound=round(100*ab[k],1),
                                                     replicates=keep,discarded=disc,
                                                     edge_RMSE=if(k==1) round(mean(g("edge_rmse")),4) else NA_real_,
                                                     edge_MedAE=if(k==1) round(median(g("edge_medae")),4) else NA_real_,
                                                     pct_runaway=if(k==1) round(100*mean(g("max_alpha")>3),1) else NA_real_,
                                                     rank_cor=if(k==1) round(mean(g("rank_cor"),na.rm=TRUE),3) else NA_real_,
                                                     rank_max=if(k==1) round(rank_ceiling(net_lead(At)),3) else NA_real_)}}
  T34<-rbindlist(S2); print(T34); fwrite(T34,file.path(OUT_DIR,"table_sim_eta.csv"))
  
  # Generates Table 5 (Study 3, robustness under a movement rule)
  S3<-list()
  for(fid in c(0.9,0.6,0.3,0.0)){
    keep<-0;det<-c();rc<-c();er<-c();nl1<-c()
    for(r in seq_len(CFG$NREP)){
      D<-simulate_follow(CFG,fid); if(min(rowSums(D$Y))<10) next
      f<-fit_sim(build_W_sim(D$mk,eta0,CFG),D$Y,CFG); if(is.null(f)) next
      keep<-keep+1; er<-c(er,mean(D$Y[,(CFG$L+1):CFG$B]))
      nl<-net_lead(f$A); det<-c(det,which.max(nl)==1); nl1<-c(nl1,nl[1])
      rc<-c(rc,suppressWarnings(cor(c(1,rep(0,CFG$M-1)),nl,method="spearman")))}
    if(!keep) next
    p<-mean(det)
    S3[[as.character(fid)]]<-data.table(fidelity=fid,replicates=keep,
                                        event_rate=round(mean(er),3),leader_detect=round(p,3),
                                        detect_mcse=round(mcse_p(p,keep),4),chance_rate=round(1/CFG$M,3),
                                        net_lead_ent1=round(mean(nl1),4),rank_cor=round(mean(rc,na.rm=TRUE),3),
                                        rank_max=round(rank_ceiling(c(1,rep(0,CFG$M-1))),3))}
  T5<-rbindlist(S3); print(T5); fwrite(T5,file.path(OUT_DIR,"table_sim_robust.csv"))
  
  # Generates Table 6 (Study 4, head-to-head against directional delay)
  S4<-list()
  for(sc in c("single","chain")){
    At<-make_A(sc,CFG$M); hk<-c();dk<-c()
    for(r in seq_len(min(CFG$NREP,100))){
      mk<-simulate_marks(CFG); W<-build_W_sim(mk,eta0,CFG)
      Y<-simulate_Y(W,At,CFG$mu,CFG,"poisson"); if(min(rowSums(Y))<10) next
      f<-fit_sim(W,Y,CFG); if(is.null(f)) next
      tl<-which(net_lead(At)>=max(net_lead(At))-1e-12)
      hk<-c(hk,which.max(net_lead(f$A))%in%tl); dk<-c(dk,which.max(delay_leader(mk,CFG))%in%tl)}
    if(!length(hk)) next
    S4[[paste0("model_",sc)]]<-data.table(dgp="Hawkes model",scenario=sc,
                                          replicates=length(hk),hawkes_correct=round(mean(hk),3),
                                          hawkes_mcse=round(mcse_p(mean(hk),length(hk)),4),
                                          delay_correct=round(mean(dk),3),delay_mcse=round(mcse_p(mean(dk),length(dk)),4),
                                          chance=round(1/CFG$M,3))}
  for(fid in c(0.9,0.6,0.3)){
    hk<-c();dk<-c()
    for(r in seq_len(min(CFG$NREP,100))){
      D<-simulate_follow(CFG,fid); if(min(rowSums(D$Y))<10) next
      f<-fit_sim(build_W_sim(D$mk,eta0,CFG),D$Y,CFG); if(is.null(f)) next
      hk<-c(hk,which.max(net_lead(f$A))==1); dk<-c(dk,which.max(delay_leader(D$mk,CFG))==1)}
    if(!length(hk)) next
    S4[[paste0("follow_",fid)]]<-data.table(dgp=paste0("following rule (fid=",fid,")"),
                                            scenario="single",replicates=length(hk),hawkes_correct=round(mean(hk),3),
                                            hawkes_mcse=round(mcse_p(mean(hk),length(hk)),4),
                                            delay_correct=round(mean(dk),3),delay_mcse=round(mcse_p(mean(dk),length(dk)),4),
                                            chance=round(1/CFG$M,3))}
  T6<-rbindlist(S4); print(T6); fwrite(T6,file.path(OUT_DIR,"table_sim_headtohead.csv"))
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Table 7: bin width and working likelihood
## --------------------------------------------------------------------
if (run("Table7")) {
  message("\n### BIN WIDTH: ONE REALIZATION, AGGREGATED ###")
  set.seed(SEED_BASE + 1010)
  
  BA <- list(
    dt0     = 0.1,                 # fine generating resolution, seconds
    B0      = 2400,                # fine bins; 240 s of flight
    factors = c(2, 4, 8, 16),      # aggregation factors -> 0.2, 0.4, 0.8, 1.6 s
    L0      = 32,                  # fine lag count; support = 3.2 s
    NREP    = 60,
    families = c("poisson", "bernoulli")
  )
  SUPPORT <- BA$L0 * BA$dt0        # 3.2 s, held fixed across resolutions
  
  agg_counts <- function(Y, f) {
    if (f == 1) return(Y)
    m <- nrow(Y); nb <- floor(ncol(Y)/f)
    out <- matrix(0, m, nb)
    for (b in seq_len(nb)) out[, b] <- rowSums(Y[, ((b-1)*f+1):(b*f), drop = FALSE])
    out
  }
  
  ## recompute the bin-level state at the coarse resolution from the SAME
  ## underlying position path
  agg_marks <- function(mk, f) {
    M <- nrow(mk$PX); nb <- floor(ncol(mk$PX)/f)
    PX <- PY <- Hd <- matrix(NA_real_, M, nb)
    for (b in seq_len(nb)) {
      idx <- ((b-1)*f+1):(b*f)
      PX[, b] <- rowMeans(mk$PX[, idx, drop = FALSE])
      PY[, b] <- rowMeans(mk$PY[, idx, drop = FALSE])
    }
    for (b in 2:nb) Hd[, b] <- atan2(PY[, b] - PY[, b-1], PX[, b] - PX[, b-1])
    Hd[, 1] <- Hd[, 2]
    list(PX = PX, PY = PY, Hd = Hd)
  }
  
  fit_at <- function(mf, Yf, cf, eta) {
    Wf <- build_W_sim(mf, eta, cf)
    A <- matrix(0, cf$M, cf$M); at <- (cf$L+1):cf$B
    for (i in seq_len(cf$M)) {
      fr <- fit_row(design_sim(Wf, Yf, i, cf)[at, , drop = FALSE], Yf[i, at], cf)
      if (is.null(fr)) return(NULL)
      A[i, ] <- fr$alpha
    }
    A
  }
  
  At   <- make_A("single", CFG$M)
  eta0 <- list(beta = CFG$beta, sigma_d = CFG$sigma_d, kappa = CFG$kappa)
  cfg0 <- modifyList(CFG, list(bin_sec = BA$dt0, B = BA$B0, L = BA$L0))
  
  ## accumulator keyed by factor and family
  key <- function(f, fam) paste(f, fam)
  ST <- list()
  for (f in BA$factors) for (fam in BA$families)
    ST[[key(f,fam)]] <- list(er=c(), rc=c(), ed=c(), rho=c(), evr=c(), disc=0)
  
  for (r in seq_len(BA$NREP)) {
    mk0 <- simulate_marks(cfg0)
    W0  <- build_W_sim(mk0, eta0, cfg0)
    ## one fine realization per family, from the same marks
    Y0 <- lapply(BA$families, function(fam) simulate_Y(W0, At, cfg0$mu, cfg0, fam))
    names(Y0) <- BA$families
    for (fam in BA$families) {
      if (min(rowSums(Y0[[fam]])) < 10) next
      for (f in BA$factors) {
        dt <- BA$dt0 * f
        L  <- max(round(SUPPORT / dt), 1)
        Yf <- agg_counts(Y0[[fam]], f)
        mf <- agg_marks(mk0, f)
        cf <- modifyList(CFG, list(bin_sec = dt, B = ncol(Yf), L = L))
        k <- key(f, fam)
        if (ncol(Yf) < L + 20) { ST[[k]]$disc <- ST[[k]]$disc + 1; next }
        A <- tryCatch(fit_at(mf, Yf, cf, eta0), error = function(e) NULL)
        if (is.null(A)) { ST[[k]]$disc <- ST[[k]]$disc + 1; next }
        m <- metrics(A, At, cf$tol_edge); off <- A; diag(off) <- 0
        ST[[k]]$er  <- c(ST[[k]]$er,  m$edge_rmse)
        ST[[k]]$rc  <- c(ST[[k]]$rc,  m$rank_cor)
        ST[[k]]$ed  <- c(ST[[k]]$ed,  sum(off > cf$tol_edge))
        ST[[k]]$rho <- c(ST[[k]]$rho, max(Mod(eigen(A, only.values = TRUE)$values)))
        ST[[k]]$evr <- c(ST[[k]]$evr, mean(Yf[, (cf$L+1):cf$B] > 0))
      }
    }
    if (r %% 10 == 0) { cat("."); flush.console() }
  }
  cat("\n")
  
  OUT <- list()
  for (f in BA$factors) for (fam in BA$families) {
    S <- ST[[key(f,fam)]]; if (!length(S$er)) next
    OUT[[key(f,fam)]] <- data.table(
      bin_sec = BA$dt0*f, L = max(round(SUPPORT/(BA$dt0*f)),1),
      family = fam, replicates = length(S$er), discarded = S$disc,
      event_rate = round(mean(S$evr),3),
      edge_RMSE = round(mean(S$er),4),
      edge_RMSE_mcse = round(mcse(S$er),5),
      median_edges = round(median(S$ed),1),
      rank_cor = round(mean(S$rc, na.rm=TRUE),3),
      rank_max = round(rank_ceiling(net_lead(At)),3),
      rho = round(mean(S$rho),3))
  }
  RES <- rbindlist(OUT); setorder(RES, bin_sec, family)
  print(RES); fwrite(RES, file.path(OUT_DIR, "binwidth_aggregated.csv"))
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Table 8: sandwich coverage under a binary response
## --------------------------------------------------------------------
if (run("Table8")) {
  message("\n### QUASI-LIKELIHOOD: SANDWICH COVERAGE ###")
  set.seed(SEED_BASE + 1020)
  eta0 <- list(beta = CFG$beta, sigma_d = CFG$sigma_d, kappa = CFG$kappa)
  SW <- list()
  for (sc in c("single","chain")) {
    At <- make_A(sc, CFG$M)
    cov_n <- cov_r <- cov_h <- c(); wid_n <- wid_r <- c(); disp <- c(); keep <- 0
    for (rr in seq_len(CFG$NREP)) {
      mk <- simulate_marks(CFG); W <- build_W_sim(mk, eta0, CFG)
      ## BERNOULLI generation: mean model correct, variance function wrong
      Y <- simulate_Y(W, At, CFG$mu, CFG, "bernoulli")
      if (min(rowSums(Y)) < 10) next
      at <- (CFG$L+1):CFG$B
      ok <- TRUE; cn <- cr <- ch <- c(); wn <- wr <- c(); dp <- c()
      for (i in seq_len(CFG$M)) {
        f <- fit_row_ql(design_sim(W,Y,i,CFG)[at,,drop=FALSE], Y[i,at], CFG)
        if (is.null(f)) { ok <- FALSE; break }
        dp <- c(dp, f$dispersion)
        for (j in which(At[i,] > 0)) {
          est <- f$alpha[j]
          for (tag in c("n","r","h")) {
            s <- switch(tag, n=f$se_naive[j], r=f$se_rob[j], h=f$se_hac[j])
            if (!is.finite(s) || s <= 0) next
            hit <- abs(est - At[i,j]) <= 1.96*s
            if (tag=="n") { cn <- c(cn,hit); wn <- c(wn, 2*1.96*s) }
            if (tag=="r") { cr <- c(cr,hit); wr <- c(wr, 2*1.96*s) }
            if (tag=="h")   ch <- c(ch,hit)
          }
        }
      }
      if (!ok) next
      keep <- keep + 1
      cov_n <- c(cov_n,cn); cov_r <- c(cov_r,cr); cov_h <- c(cov_h,ch)
      wid_n <- c(wid_n,wn); wid_r <- c(wid_r,wr); disp <- c(disp,dp)
    }
    if (!keep) next
    SW[[sc]] <- data.table(
      scenario = sc, replicates = keep,
      mean_dispersion = round(mean(disp, na.rm=TRUE), 3),
      cov_naive = round(mean(cov_n, na.rm=TRUE), 3),
      mcse_naive = round(mcse_p(mean(cov_n,na.rm=TRUE), sum(is.finite(cov_n))), 4),
      cov_robust = round(mean(cov_r, na.rm=TRUE), 3),
      mcse_robust = round(mcse_p(mean(cov_r,na.rm=TRUE), sum(is.finite(cov_r))), 4),
      cov_hac = round(mean(cov_h, na.rm=TRUE), 3),
      width_ratio = round(mean(wid_r,na.rm=TRUE)/mean(wid_n,na.rm=TRUE), 3))
  }
  TSW <- rbindlist(SW); print(TSW)
  fwrite(TSW, file.path(OUT_DIR,"table_sandwich_coverage.csv"))
  cat("\nEXPECTED PATTERN: with Bernoulli generation Var(Y)=lambda(1-lambda) < lambda,\n")
  cat("so the naive Poisson interval is too WIDE and OVER-covers. The robust\n")
  cat("interval should sit at nominal 0.95 with a width ratio below one. If the\n")
  cat("naive coverage is at or below nominal, check the generation family.\n")
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: rho(R), intensity shares, dAIC for free diagonal

## --------------------------------------------------------------------
##  Section 6: threshold-free edge selection
## --------------------------------------------------------------------
if (run("Sec6Edges")) {
  message("\n### EDGE SELECTION WITHOUT AN ORACLE THRESHOLD ###")
  set.seed(SEED_BASE + 1030)
  
  ES <- list(
    NREP      = 200,
    R_BOOT    = 200,   # bootstrap replicates per replicate, for the interval rule
    N_SUB     = 100,   # subsamples for stability selection
    SUB_FRAC  = 0.5,   # fraction of bins per subsample
    PI_THRESH = 0.60,  # selection frequency required by the stability rule
    REL_FRAC  = 0.10   # fraction of the largest off-diagonal for the relative rule
  )
  
  ## --- threshold-free summaries ----------------------------------------
  ## AUC by the rank-sum identity; average precision by sweeping cutoffs.
  auc_of <- function(score, lab) {
    n1 <- sum(lab); n0 <- sum(!lab)
    if (n1 == 0 || n0 == 0) return(NA_real_)
    r <- rank(score, ties.method = "average")
    (sum(r[lab]) - n1*(n1+1)/2) / (n1*n0)
  }
  avg_prec <- function(score, lab) {
    n1 <- sum(lab); if (n1 == 0) return(NA_real_)
    o <- order(score, decreasing = TRUE); lab <- lab[o]
    tp <- cumsum(lab); prec <- tp/seq_along(lab)
    sum(prec[lab])/n1
  }
  
  ## --- one replicate ----------------------------------------------------
  one_rep <- function(At, cfg, eta) {
    mk <- simulate_marks(cfg); W <- build_W_sim(mk, eta, cfg)
    Y  <- simulate_Y(W, At, cfg$mu, cfg, "poisson")
    if (min(rowSums(Y)) < 10) return(NULL)
    at <- (cfg$L+1):cfg$B
    Ahat <- fit_sim(W, Y, cfg); if (is.null(Ahat)) return(NULL)
    A <- Ahat$A
    
    ## bootstrap intervals: resample bins with replacement within each target row
    LO <- matrix(NA_real_, cfg$M, cfg$M)
    boot <- array(NA_real_, c(ES$R_BOOT, cfg$M, cfg$M))
    Sall <- lapply(seq_len(cfg$M), function(i) design_sim(W, Y, i, cfg))
    for (r in seq_len(ES$R_BOOT)) {
      idx <- sample(at, length(at), replace = TRUE)
      for (i in seq_len(cfg$M)) {
        fr <- fit_row(Sall[[i]][idx, , drop = FALSE], Y[i, idx], cfg)
        if (!is.null(fr)) boot[r, i, ] <- fr$alpha
      }
    }
    for (i in seq_len(cfg$M)) for (j in seq_len(cfg$M))
      LO[i,j] <- quantile(boot[, i, j], 0.025, na.rm = TRUE, names = FALSE)
    
    ## stability selection: refit on subsamples of the bins
    FREQ <- matrix(0, cfg$M, cfg$M)
    nsub <- floor(ES$SUB_FRAC*length(at))
    for (s in seq_len(ES$N_SUB)) {
      idx <- sort(sample(at, nsub))
      for (i in seq_len(cfg$M)) {
        fr <- fit_row(Sall[[i]][match(idx, at), , drop = FALSE], Y[i, idx], cfg)
        if (!is.null(fr)) FREQ[i, ] <- FREQ[i, ] + (fr$alpha > cfg$tol_edge)
      }
    }
    FREQ <- FREQ/ES$N_SUB
    
    ## assemble off-diagonal vectors
    keep <- which(row(At) != col(At))
    list(score = A[keep], truth = At[keep] > 0,
         lo = LO[keep], freq = FREQ[keep],
         relmax = if (max(A[keep]) > 0) A[keep]/max(A[keep]) else A[keep]*0,
         oracle_thr = 0.5*min(At[At > 0 & row(At) != col(At)]))
  }
  
  prf <- function(sel, truth) {
    TP <- sum(sel & truth); FP <- sum(sel & !truth); FN <- sum(!sel & truth)
    c(sens = if (TP+FN > 0) TP/(TP+FN) else NA_real_,
      prec = if (TP+FP > 0) TP/(TP+FP) else NA_real_)
  }
  
  OUT <- list()
  for (sc in c("single","chain","coleader","subgroup")) {
    At  <- make_A(sc, CFG$M)
    eta <- list(beta = CFG$beta, sigma_d = CFG$sigma_d, kappa = CFG$kappa)
    AUC <- AP <- c(); RUL <- list()
    for (r in seq_len(ES$NREP)) {
      z <- tryCatch(one_rep(At, CFG, eta), error = function(e) NULL)
      if (is.null(z)) next
      AUC <- c(AUC, auc_of(z$score, z$truth))
      AP  <- c(AP,  avg_prec(z$score, z$truth))
      RUL[[length(RUL)+1L]] <- rbind(
        oracle    = prf(z$score  >  z$oracle_thr,  z$truth),
        bootstrap = prf(z$lo     >  0,             z$truth),
        stability = prf(z$freq   >= ES$PI_THRESH,  z$truth),
        relative  = prf(z$relmax >= ES$REL_FRAC,   z$truth))
      if (r %% 25 == 0) { cat("."); flush.console() }
    }
    if (!length(RUL)) next
    for (rule in c("oracle","bootstrap","stability","relative")) {
      se <- vapply(RUL, function(m) m[rule,"sens"], numeric(1))
      pr <- vapply(RUL, function(m) m[rule,"prec"], numeric(1))
      OUT[[paste(sc,rule)]] <- data.table(
        scenario = sc, rule = rule, replicates = length(RUL),
        sensitivity = round(mean(se, na.rm=TRUE), 3),
        sens_mcse   = round(mcse(se), 4),
        precision   = round(mean(pr, na.rm=TRUE), 3),
        prec_mcse   = round(mcse(pr), 4),
        AUC = if (rule=="oracle") round(mean(AUC, na.rm=TRUE), 3) else NA_real_,
        avg_precision = if (rule=="oracle") round(mean(AP, na.rm=TRUE), 3) else NA_real_)
    }
    cat(sprintf("\n  %-9s AUC %.3f  AP %.3f\n", sc, mean(AUC,na.rm=TRUE), mean(AP,na.rm=TRUE)))
  }
  RES <- rbindlist(OUT); setorder(RES, scenario, rule)
  print(RES); fwrite(RES, file.path(OUT_DIR, "edge_selection.csv"))
  
  message("  Complete")
}
message("\nDone. Outputs written to ", normalizePath(OUT_DIR))
sessionInfo()

## --------------------------------------------------------------------
##  Section 6: shared group-level driver
## --------------------------------------------------------------------
if (run("Sec6Driver")) {
  message("\n### SHARED DRIVER SIMULATION ###")
  set.seed(SEED_BASE + 1040)
  
  
  ARMS <- list(
    "driver only"       = list(drive=TRUE,  coupled=FALSE),
    "driver + coupling" = list(drive=TRUE,  coupled=TRUE ),
    "coupling only"     = list(drive=FALSE, coupled=TRUE ))
  
  OUT <- list()
  for (arm in names(ARMS)) for (use_gamma in c(FALSE, TRUE)) {
    cfgA <- CFG
    At <- if (ARMS[[arm]]$coupled) make_A("single", cfgA$M) else
    { Z <- matrix(0, cfgA$M, cfgA$M); diag(Z) <- 0.40; Z }
    eta0 <- list(beta=cfgA$beta, sigma_d=cfgA$sigma_d, kappa=cfgA$kappa)
    ro <- ri <- c(); no_p <- no_c <- ni_s <- list(); keep <- 0
    for (r in seq_len(SD$NREP)) {
      mk <- simulate_marks(cfgA); W <- build_W_sim(mk, eta0, cfgA)
      drv <- if (ARMS[[arm]]$drive) make_driver(cfgA$B, cfgA) else rep(1, cfgA$B)
      Y <- simulate_arm(W, At, drv, cfgA)
      if (min(rowSums(Y)) < 10) next
      f0 <- fit_arm(W, Y, cfgA, use_gamma); if (is.null(f0)) next
      keep <- keep + 1
      ## test bird 1, the designated leader in the coupled arms
      b <- 1
      ro <- c(ro, f0$out[b]); ri <- c(ri, f0$inc[b])
      pn <- cn <- numeric(SD$N_PERM)
      for (k in seq_len(SD$N_PERM)) {
        Yp <- Y; Yp[b,] <- perm_row(Y[b,])
        fp <- fit_arm(W, Y, cfgA, use_gamma, Yover=Yp)
        pn[k] <- if (is.null(fp)) NA_real_ else fp$out[b]
        cv <- circ_row(Y[b,], cfgA$L)
        if (is.null(cv)) { cn[k] <- NA_real_ } else {
          Yc <- Y; Yc[b,] <- cv
          fc <- fit_arm(W, Y, cfgA, use_gamma, Yover=Yc)
          cn[k] <- if (is.null(fc)) NA_real_ else fc$out[b]
        }
      }
      no_p[[keep]] <- pn[is.finite(pn)]; no_c[[keep]] <- cn[is.finite(cn)]
    }
    if (!keep) next
    pp <- mapply(emp_p, ro, no_p); pc <- mapply(emp_p, ro, no_c)
    nmp <- vapply(no_p, function(x) if(length(x)) median(x) else NA_real_, numeric(1))
    nmc <- vapply(no_c, function(x) if(length(x)) median(x) else NA_real_, numeric(1))
    OUT[[paste(arm, use_gamma)]] <- data.table(
      arm = arm, gamma = if (use_gamma) "free" else "fixed at 0",
      replicates = keep,
      real_out = round(median(ro), 4),
      null_perm = round(median(nmp, na.rm=TRUE), 4),
      null_circ = round(median(nmc, na.rm=TRUE), 4),
      pct_above_perm = round(100*mean(ro > nmp, na.rm=TRUE), 1),
      pct_p_lt_05 = round(100*mean(pp < 0.05, na.rm=TRUE), 1),
      floor_p = round(1/(SD$N_PERM+1), 3),
      fisher_p = signif(fisher(pp), 3))
    cat(sprintf("  %-18s gamma %-10s  real %.4f  null %.4f  %%above %.0f\n",
                arm, if (use_gamma) "free" else "fixed", median(ro),
                median(nmp, na.rm=TRUE), 100*mean(ro > nmp, na.rm=TRUE)))
  }
  RES <- rbindlist(OUT)
  print(RES); fwrite(RES, file.path(OUT_DIR, "shared_driver_simulation.csv"))
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Section 6: calibrated coupling sweep
## --------------------------------------------------------------------
if (run("Sec6Calib")) {
  message("\n### CALIBRATED COUPLING SWEEP ###")
  set.seed(SEED_BASE + 1045)
  
  ALPHAS    <- c(0.00, 0.03, 0.064, 0.10, 0.15, 0.30)
  PIG_RATIO <- 3.4     # applied value, real/null under the group-activity model
  
  CAL <- list()
  for (ax in ALPHAS) for (use_gamma in c(FALSE, TRUE)) {
    cfgA <- CFG
    At <- matrix(0, cfgA$M, cfgA$M); diag(At) <- 0.40
    if (ax > 0) At[2:cfgA$M, 1] <- ax           # single leader at strength ax
    eta0 <- list(beta = cfgA$beta, sigma_d = cfgA$sigma_d, kappa = cfgA$kappa)
    ro <- c(); nl <- list(); keep <- 0
    for (r in seq_len(SD$NREP)) {
      mk  <- simulate_marks(cfgA)
      W   <- build_W_sim(mk, eta0, cfgA)
      drv <- make_driver(cfgA$B, cfgA)          # driver present in every cell
      Y   <- simulate_arm(W, At, drv, cfgA)
      if (min(rowSums(Y)) < 10) next
      f0 <- fit_arm(W, Y, cfgA, use_gamma); if (is.null(f0)) next
      keep <- keep + 1; b <- 1
      ro <- c(ro, f0$out[b])
      pn <- numeric(SD$N_PERM)
      for (k in seq_len(SD$N_PERM)) {
        Yp <- Y; Yp[b, ] <- perm_row(Y[b, ])
        fp <- fit_arm(W, Y, cfgA, use_gamma, Yover = Yp)
        pn[k] <- if (is.null(fp)) NA_real_ else fp$out[b]
      }
      nl[[keep]] <- pn[is.finite(pn)]
    }
    if (!keep) next
    nm <- vapply(nl, function(x) if (length(x)) median(x) else NA_real_, numeric(1))
    pp <- mapply(emp_p, ro, nl)
    rr <- median(ro, na.rm = TRUE) / median(nm, na.rm = TRUE)
    CAL[[paste(ax, use_gamma)]] <- data.table(
      alpha = ax, gamma = if (use_gamma) "free" else "fixed at 0",
      replicates = keep,
      real_out  = round(median(ro, na.rm = TRUE), 4),
      null_med  = round(median(nm, na.rm = TRUE), 4),
      ratio     = round(rr, 2),
      pct_above = round(100*mean(ro > nm, na.rm = TRUE), 1),
      fisher_p  = signif(fisher(pp), 3))
    cat(sprintf("  alpha %.3f  gamma %-10s  ratio %5.2f  %%above %5.1f\n",
                ax, if (use_gamma) "free" else "fixed", rr,
                100*mean(ro > nm, na.rm = TRUE)))
  }
  RES <- rbindlist(CAL)
  print(RES); fwrite(RES, file.path(OUT_DIR, "calibrated_coupling.csv"))
  
  gf <- RES[gamma == "free"]
  if (nrow(gf) > 1) {
    cat(sprintf("\n=== applied ratio %.2f against the sweep (gamma free) ===\n",
                PIG_RATIO))
    for (i in seq_len(nrow(gf))) {
      mark <- if (abs(gf$ratio[i] - PIG_RATIO) ==
                  min(abs(gf$ratio - PIG_RATIO))) "   <-- closest" else ""
      cat(sprintf("  alpha %.3f -> ratio %5.2f%s\n", gf$alpha[i], gf$ratio[i], mark))
    }
    br <- gf$ratio[gf$alpha == 0]
    cat(sprintf("\n  coupling-free ratio: %.2f\n", br))
    cat(sprintf("  median non-zero off-diagonal estimated in the application: 0.064\n"))
    cat("\n  Report the ratio at alpha = 0 and at alpha = 0.064 in the text, and\n")
    cat("  state where the applied value falls between them. Do NOT report the\n")
    cat("  proportion above the null from this sweep alongside the applied\n")
    cat("  figure: the simulated flocks are homogeneous while the releases are\n")
    cat("  not, so that quantity is not comparable between the two.\n")
  }
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Section 7: intensity decomposition and rho(R)
## --------------------------------------------------------------------
if (run("Sec7Fits")) {
  message("\n### PIGEON REAL FITS ###")
  set.seed(SEED_BASE + 1050)
  P<-prep_all()
  FREE<-FIX<-list()
  for(fk in names(P)){
    M<-P[[fk]]; DZ<-build_W(M,CFG)
    FREE[[fk]]<-tryCatch(fit_flight(DZ,M,CFG,want_se=TRUE,want_Ew=TRUE),error=function(e)NULL)
    ## alpha_ii = 0 comparison: zero the diagonal column from each design
    DZ0<-DZ
    f0<-tryCatch({al<-matrix(0,M$m,M$m); ll<-0; np<-0; sh<-matrix(NA_real_,M$m,3)
    for(i in seq_len(M$m)){at<-DZ$rows[[i]]; if(length(at)<2) next
    y<-M$Y[i,at]; if(sum(y)<CFG$min_events) next
    S<-design_i(DZ,M$Y,i,CFG)[at,,drop=FALSE]; S[,i]<-0
    r<-fit_row(S,y,CFG); if(is.null(r)) next
    al[i,]<-r$alpha; ll<-ll+r$loglik; np<-np+r$npar
    con<-sweep(S,2,r$alpha,`*`); tot<-mean(r$mu+rowSums(con))
    if(is.finite(tot)&&tot>0) sh[i,]<-c(r$mu/tot,0,mean(rowSums(con))/tot)}
    list(loglik=ll,npar=np,shares=colMeans(sh,na.rm=TRUE))},error=function(e)NULL)
    FIX[[fk]]<-f0}
  FREE<-Filter(Negate(is.null),FREE); FIX<-Filter(Negate(is.null),FIX)
  fv<-function(L,k) vapply(L,`[[`,numeric(1),k)
  shF<-rowMeans(vapply(FREE,`[[`,numeric(3),"shares"),na.rm=TRUE)
  shX<-rowMeans(vapply(FIX ,`[[`,numeric(3),"shares"),na.rm=TRUE)
  rhoR<-fv(FREE,"rho_R")
  RF<-data.table(quantity=c("median rho(R)","max rho(R)","rho(R)>=1 count",
                            "rho(R) 95% CI lower","rho(R) 95% CI upper","median rho(A) off-diagonal",
                            "median alpha_ii","median edges/flight","background % (alpha_ii free)",
                            "self % (alpha_ii free)","cross % (alpha_ii free)",
                            "background % (alpha_ii=0)","cross % (alpha_ii=0)","dAIC free vs fixed diagonal"),
                 value=c(round(median(rhoR),3),round(max(rhoR),3),sum(rhoR>=1),
                         round(quantile(rhoR,.025,names=FALSE),3),round(quantile(rhoR,.975,names=FALSE),3),
                         round(median(fv(FREE,"rho_off")),4),
                         round(median(unlist(lapply(FREE,`[[`,"self")),na.rm=TRUE),3),
                         round(median(fv(FREE,"n_edges")),0),
                         round(100*shF[1],1),round(100*shF[2],1),round(100*shF[3],1),
                         round(100*shX[1],1),round(100*shX[3],1),
                         round((2*sum(fv(FIX,"npar"))-2*sum(fv(FIX,"loglik")))-
                                 (2*sum(fv(FREE,"npar"))-2*sum(fv(FREE,"loglik"))),1)))
  print(RF); fwrite(RF,file.path(OUT_DIR,"pigeon_real_fits.csv"))
  saveRDS(FREE,file.path(OUT_DIR,"ckpt","free_fits.rds"))
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Table 9 (null battery) and Section 8.3 p-values

## --------------------------------------------------------------------
##  Section 7: cluster bootstrap intervals
## --------------------------------------------------------------------
if (run("Sec7Boot")) {
  message("\n### CLUSTER BOOTSTRAP OVER FLOCK-RELEASES ###")
  set.seed(SEED_BASE + 1060)
  P <- prep_all()
  
  ## fit every release once; the bootstrap then resamples these fitted objects
  FITS <- list()
  for (fk in names(P)) {
    G <- P[[fk]]
    f <- tryCatch(fit_flight(build_W(G,CFG), G, CFG, want_Ew = TRUE), error=function(e) NULL)
    if (!is.null(f)) FITS[[fk]] <- f
  }
  cat("fitted", length(FITS), "releases\n")
  
  ## summaries computed on a COLLECTION of releases; these are the estimands
  summarise_releases <- function(FL) {
    if (!length(FL)) return(NULL)
    off  <- unlist(lapply(FL, function(z) z$off[z$off > CFG$tol_edge]))
    self <- unlist(lapply(FL, `[[`, "self"))
    out  <- unlist(lapply(FL, `[[`, "out"))
    inc  <- unlist(lapply(FL, `[[`, "inc"))
    sh   <- rowMeans(vapply(FL, `[[`, numeric(3), "shares"), na.rm = TRUE)
    rhoR <- vapply(FL, `[[`, numeric(1), "rho_R")
    c(median_off_alpha = median(off, na.rm=TRUE),
      median_self      = median(self, na.rm=TRUE),
      median_out       = median(out,  na.rm=TRUE),
      median_inc       = median(inc,  na.rm=TRUE),
      share_mu         = 100*sh[1],
      share_self       = 100*sh[2],
      share_cross      = 100*sh[3],
      median_rho_R     = median(rhoR, na.rm=TRUE),
      max_rho_R        = max(rhoR, na.rm=TRUE),
      median_edges     = median(vapply(FL,`[[`,numeric(1),"n_edges")),
      frac_rho_gt_95 = mean(rhoR > 0.95),
      frac_rho_gt_99 = mean(rhoR > 0.99),
      q95_rho_R      = quantile(rhoR, 0.95, na.rm=TRUE))
  }
  hat <- summarise_releases(FITS)
  nm  <- names(hat)
  R   <- 2000                      # cheap: no refitting, only resampling releases
  BS  <- matrix(NA_real_, R, length(hat), dimnames = list(NULL, nm))
  keys <- names(FITS)
  for (r in seq_len(R)) {
    idx <- sample(keys, length(keys), replace = TRUE)
    BS[r,] <- summarise_releases(FITS[idx])
  }
  CB <- data.table(
    quantity = nm,
    estimate = round(hat, 4),
    lower    = round(apply(BS, 2, quantile, .025, na.rm=TRUE), 4),
    upper    = round(apply(BS, 2, quantile, .975, na.rm=TRUE), 4),
    boot_se  = round(apply(BS, 2, sd, na.rm=TRUE), 4))
  print(CB); fwrite(CB, file.path(OUT_DIR,"cluster_bootstrap.csv"))
  ## Sensitivity: do the near-critical releases drive the decomposition?
  ## 26 of 116 releases have rho(R-hat) > 0.95. If the intensity shares are
  ## unchanged when those are excluded, the decomposition does not depend on them.
  rhoR_all <- vapply(FITS, `[[`, numeric(1), "rho_R")
  FITS_sub <- FITS[rhoR_all <= 0.95]
  cat(sprintf("\nexcluding rho(R) > 0.95: %d of %d releases retained\n",
              length(FITS_sub), length(FITS)))
  sub_hat <- summarise_releases(FITS_sub)
  
  ## bootstrap the retained subset the same way, so the comparison is like for like
  BS2 <- matrix(NA_real_, R, length(sub_hat), dimnames = list(NULL, names(sub_hat)))
  k2 <- names(FITS_sub)
  for (r in seq_len(R)) BS2[r,] <- summarise_releases(FITS_sub[sample(k2, length(k2), replace = TRUE)])
  CB2 <- data.table(
    quantity = names(sub_hat),
    full     = round(hat, 4),
    excl_hi  = round(sub_hat, 4),
    lower    = round(apply(BS2, 2, quantile, .025, na.rm = TRUE), 4),
    upper    = round(apply(BS2, 2, quantile, .975, na.rm = TRUE), 4))
  print(CB2); fwrite(CB2, file.path(OUT_DIR, "cluster_bootstrap_subcritical.csv"))
  cat("\nNOTE ON WHAT IS AND IS NOT BOOTSTRAPPED HERE.\n")
  cat("  Releases are resampled; the fit within a release is not recomputed.\n")
  cat("  This propagates variation BETWEEN releases, which is the level at which\n")
  cat("  releases are independent, and is the correct unit given that birds recur\n")
  cat("  across releases and are dependent within them. It does NOT propagate\n")
  cat("  within-release estimation error. If both are wanted, refit inside the\n")
  cat("  loop; that is roughly R times the cost of a single pass and is only\n")
  cat("  worth it if the within-release component is material, which can be\n")
  cat("  checked by comparing boot_se above with the per-release spread.\n")
  message("  Complete")
}
## ---------------------------------------------------------------------
#  Generates Section 8.4: event-definition sensitivity (thresholds and types)

## --------------------------------------------------------------------
##  Section 7: residual diagnostics
## --------------------------------------------------------------------
if (run("Sec7Resid")) {
  message("\n### RESIDUAL DIAGNOSTICS ###")
  set.seed(SEED_BASE + 1070)
  P<-prep_all(); FL<-list()
  for(fk in names(P)){G<-P[[fk]]
  f<-tryCatch(fit_flight(build_W(G,CFG),G,CFG,want_fit=TRUE),error=function(e)NULL)
  if(is.null(f)||is.null(f$fitted)) next
  r<-f$fitted; r[,frkey:=fk]; FL[[fk]]<-r}
  D<-rbindlist(FL)
  D[,r_pois:=qnorm(pmin(pmax(runif(.N,ppois(y-1,lambda),ppois(y,lambda)),1e-10),1-1e-10))]
  D[,p_:=pmin(pmax(lambda,1e-10),1-1e-10)]
  D[,r_bern:=qnorm(pmin(pmax(runif(.N,ifelse(y==0,0,1-p_),ifelse(y==0,1-p_,1)),1e-10),1-1e-10))]
  D[,pear:=(y-lambda)/sqrt(pmax(lambda,1e-12))]
  np<-sum(vapply(P,function(g) g$m,numeric(1)))
  DG<-data.table(quantity=c("RQR Poisson","RQR Bernoulli"),
                 mean=round(c(mean(D$r_pois),mean(D$r_bern)),4),
                 sd=round(c(sd(D$r_pois),sd(D$r_bern)),4))
  print(DG)
  cat(sprintf("dispersion %.3f | obs zeros %.4f | Poisson-pred %.4f | Bernoulli-pred %.4f\n",
              sum(D$pear^2)/(nrow(D)-np),mean(D$y==0),mean(exp(-D$lambda)),mean(1-pmin(D$lambda,1))))
  ACF<-D[,{r<-r_bern[order(bin)]
  if(length(r)<CFG$L+5) .(lag=integer(),rho=numeric())
  else .(lag=1:CFG$L,rho=vapply(1:CFG$L,function(k)
    suppressWarnings(cor(r[-seq_len(k)],r[seq_len(length(r)-k)])),numeric(1)))},
  by=.(frkey,bird)][is.finite(rho),.(mean_rho=round(mean(rho),4),n=.N),by=lag][order(lag)]
  print(ACF)
  D[,dec:=cut(lambda,breaks=quantile(lambda,seq(0,1,length.out=CFG$n_decile+1)),
              include.lowest=TRUE,labels=FALSE)]
  CAL<-D[,.(n=.N,fitted=round(mean(lambda),4),observed=round(mean(y),4)),by=dec][order(dec)]
  CAL[,diff:=round(observed-fitted,4)]; print(CAL)
  D[,tdec:=cut(t_rel,breaks=seq(0,1,length.out=11),include.lowest=TRUE,labels=FALSE)]
  BT<-D[,.(n=.N,mean_r=round(mean(r_bern),4)),by=tdec][order(tdec)]
  pdf(file.path(OUT_DIR,"figure_diagnostics.pdf"),width=9,height=8)
  par(mfrow=c(2,2),mar=c(4.2,4.2,2.5,1))
  s<-D[sample(.N,min(20000,.N))]
  qqnorm(s$r_bern,pch=16,cex=.25,col=rgb(0,0,0,.3),main="(a) Normal Q-Q"); qqline(s$r_bern,col=2)
  plot(CAL$dec,D[,.(m=mean(r_bern)),by=dec][order(dec)]$m,type="b",pch=16,ylim=c(-.5,.5),
       xlab="Decile of fitted mean",ylab="Mean residual",main="(b) Residual vs fitted"); abline(h=0,lty=2)
  plot(BT$tdec,BT$mean_r,type="b",pch=16,ylim=c(-.5,.5),
       xlab="Decile of time in flight",ylab="Mean residual",main="(c) Residual vs time"); abline(h=0,lty=2)
  plot(ACF$lag,ACF$mean_rho,type="h",lwd=3,xlab="Lag (bins)",ylab="Mean autocorrelation",
       main="(d) Residual ACF"); abline(h=0); abline(h=c(-2,2)/sqrt(500),lty=2,col="grey40")
  dev.off()
  fwrite(DG,file.path(OUT_DIR,"diagnostics_summary.csv"))
  fwrite(ACF,file.path(OUT_DIR,"diagnostics_acf.csv"))
  fwrite(CAL,file.path(OUT_DIR,"diagnostics_calibration.csv"))
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: cross-bird check on the 2 s residual oscillation

## --------------------------------------------------------------------
##  Section 7: 2 s oscillation, within and cross bird
## --------------------------------------------------------------------
if (run("Sec7Lag")) {
  message("\n### 2 s OSCILLATION: WITHIN- AND CROSS-BIRD ###")
  set.seed(SEED_BASE + 1080)
  P<-prep_all()
  lagcor<-function(a,b,k){n<-length(a); if(k>=n) return(NA_real_)
  x<-a[(k+1):n];y<-b[1:(n-k)]; ok<-is.finite(x)&is.finite(y)
  if(sum(ok)<50) return(NA_real_); suppressWarnings(cor(x[ok],y[ok]))}
  cshift<-function(v,k){n<-length(v); v[((seq_len(n)-1+k)%%n)+1]}
  W_<-list();X_<-list();XP<-list()
  for(fk in names(P)){
    G<-P[[fk]]; f<-tryCatch(fit_flight(build_W(G,CFG),G,CFG,want_fit=TRUE),error=function(e)NULL)
    if(is.null(f)||is.null(f$fitted)) next
    Rm<-matrix(NA_real_,G$m,G$B)
    ft<-f$fitted; ft[,p_:=pmin(pmax(lambda,1e-10),1-1e-10)]
    ft[,r:=qnorm(pmin(pmax(runif(.N,ifelse(y==0,0,1-p_),ifelse(y==0,1-p_,1)),1e-10),1-1e-10))]
    for(k in seq_len(nrow(ft))) Rm[ft$bird[k],ft$bin[k]]<-ft$r[k]
    ok<-which(rowSums(is.finite(Rm))>100); if(length(ok)<2) next
    for(i in ok) for(k in 1:12){r<-lagcor(Rm[i,],Rm[i,],k)
    if(is.finite(r)) W_[[length(W_)+1L]]<-data.table(lag=k,rho=r)}
    for(i in ok) for(j in ok){ if(i==j) next
      for(k in c(0,1,4,5,6,10)){r<-lagcor(Rm[i,],Rm[j,],k)
      if(is.finite(r)) X_[[length(X_)+1L]]<-data.table(lag=k,rho=r)}}
    for(p in 1:20){Rp<-Rm; for(i in ok) Rp[i,]<-cshift(Rm[i,],sample(ncol(Rm),1))
    for(i in ok) for(j in ok){ if(i==j) next
      for(k in c(0,5)){r<-lagcor(Rp[i,],Rp[j,],k)
      if(is.finite(r)) XP[[length(XP)+1L]]<-data.table(lag=k,rho=r)}}}}
  Wd<-rbindlist(W_);Xd<-rbindlist(X_);XPd<-rbindlist(XP)
  cat("within-bird ACF:\n"); print(Wd[,.(mean_rho=round(mean(rho),4),n=.N),by=lag][order(lag)])
  cat("cross-bird, real:\n"); print(Xd[,.(mean_rho=round(mean(rho),4),n=.N),by=lag][order(lag)])
  cat("cross-bird, circular-shift control:\n")
  print(XPd[,.(mean_rho=round(mean(rho),4),n=.N),by=lag][order(lag)])
  fwrite(Wd[,.(mean_rho=mean(rho)),by=lag],file.path(OUT_DIR,"lag5_within.csv"))
  fwrite(Xd[,.(mean_rho=mean(rho)),by=lag],file.path(OUT_DIR,"lag5_cross.csv"))
  fwrite(XPd[,.(mean_rho=mean(rho)),by=lag],file.path(OUT_DIR,"lag5_control.csv"))
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: origin tests for the 2 s oscillation

## --------------------------------------------------------------------
##  Section 7: raw-signal spectrum before binning
## --------------------------------------------------------------------
if (run("Sec7Origin")) {
  message("\n### 2 s OSCILLATION: ORIGIN ###")
  set.seed(SEED_BASE + 1090)
  need_pigeons()
  A<-copy(PIG); A[,dt:=as.numeric(time-shift(time)),by=.(frkey,bird)]
  dts<-A[is.finite(dt)&dt>0&dt<5,dt]
  cat(sprintf("median inter-fix %.4f s | %% within 10%% of nominal %.1f | %% dropped %.2f\n",
              median(dts),100*mean(abs(dts-1/CFG$fs)<0.1/CFG$fs),100*mean(dts>2/CFG$fs)))
  A[,idx:=seq_len(.N),by=.(frkey,bird)]; A[,gap:=as.integer(is.finite(dt)&dt>1.5/CFG$fs)]
  print(A[,.(p=round(mean(gap,na.rm=TRUE),5)),by=.(mod10=idx%%10)][order(mod10)])
  spec_at<-function(v,fs,period){v<-v[is.finite(v)];n<-length(v)
  if(n<256) return(NULL); v<-v-mean(v)
  sp<-spec.pgram(ts(v,frequency=fs),taper=0.1,detrend=TRUE,fast=TRUE,plot=FALSE)
  f0<-1/period; nr<-which(abs(sp$freq-f0)<0.05)
  bs<-which(abs(sp$freq-f0)>=0.05&abs(sp$freq-f0)<0.5)
  if(!length(nr)||!length(bs)) return(NULL)
  data.table(peak_ratio=max(sp$spec[nr])/median(sp$spec[bs]),
             peak_freq=sp$freq[nr][which.max(sp$spec[nr])])}
  set.seed(1); fl<-sample(unique(A$frkey),min(60,uniqueN(A$frkey))); SP<-list()
  for(fk in fl){d<-A[frkey==fk&speed>=CFG$speed_kmh_fly]
  for(b in unique(d$bird)){bb<-d[bird==b]; if(nrow(bb)<300) next
  dd<-sqrt(diff(bb$x)^2+diff(bb$y)^2)
  s1<-spec_at(bb$speed,CFG$fs,2.0); s2<-spec_at(dd,CFG$fs,2.0)
  if(!is.null(s1)) SP[[length(SP)+1L]]<-data.table(sig="raw speed",s1)
  if(!is.null(s2)) SP[[length(SP)+1L]]<-data.table(sig="raw displacement",s2)}}
  SPd<-rbindlist(SP)
  print(SPd[,.(n=.N,median_peak_ratio=round(median(peak_ratio),3),
               median_peak_freq=round(median(peak_freq),4)),by=sig])
  fwrite(SPd,file.path(OUT_DIR,"origin_spectrum.csv"))
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 6 note: predictable vs contemporaneous target marks

## --------------------------------------------------------------------
##  Table 9: null battery
## --------------------------------------------------------------------
if (run("Table9")) {
  message("\n### TABLE 9: NULL BATTERY ###")
  set.seed(SEED_BASE + 1100)
  P<-prep_all()
  bird_fl<-list()
  for(fk in names(P)) for(b in P[[fk]]$birds){
    k<-paste(flock_of(fk),b,sep="|"); bird_fl[[k]]<-c(bird_fl[[k]],fk)}
  for(fk in names(P)){
    ck<-file.path(OUT_DIR,"ckpt",paste0("t9_",gsub("[^A-Za-z0-9]","_",fk),".rds"))
    if(file.exists(ck)) next
    G<-P[[fk]]; fl<-flock_of(fk); DZ<-build_W(G,CFG); f0<-fit_flight(DZ,G,CFG)
    rows<-list()
    for(bi in seq_along(G$birds)){
      b<-G$birds[bi]
      dn<-setdiff(bird_fl[[paste(fl,b,sep="|")]],fk); inc<-numeric(0)
      if(length(dn)){ if(length(dn)>CFG$N_SWAP) dn<-sample(dn,CFG$N_SWAP)
      for(d in dn){
        dr<-tryCatch(donor_reg(G,P[[d]]$reg,b,CFG),error=function(e)NULL)
        if(is.null(dr)) next
        rs<-rbind(G$reg[bird!=b,.(bird,bin,x,y,disp,heading,event)],dr)
        Ms<-mats_from_reg(rs,G$birds,G$B)
        fs<-tryCatch(fit_flight(build_W(Ms,CFG),Ms,CFG),error=function(e)NULL)
        if(!is.null(fs)) inc<-c(inc,fs$inc[bi])}}
      pn<-cn<-mn<-numeric(CFG$N_PERM)
      for(r in seq_len(CFG$N_PERM)){
        pn[r]<-out_null(DZ,G,bi,"perm",CFG)
        cn[r]<-out_null(DZ,G,bi,"circ",CFG)
        mn[r]<-out_null(DZ,G,bi,"permmark",CFG)}
      rows[[length(rows)+1L]]<-list(fk=fk,bird=b,real_inc=f0$inc[bi],real_out=f0$out[bi],
                                    inc_swap=inc,out_perm=pn[is.finite(pn)],out_circ=cn[is.finite(cn)],
                                    out_pm=mn[is.finite(mn)])}
    saveRDS(rows,ck); cat(".")}
  cat("\n")
  rows<-unlist(lapply(list.files(file.path(OUT_DIR,"ckpt"),pattern="^t9_",full.names=TRUE),
                      readRDS),recursive=FALSE)
  ri<-vapply(rows,`[[`,numeric(1),"real_inc"); ro<-vapply(rows,`[[`,numeric(1),"real_out")
  pool<-function(f) unlist(lapply(rows,`[[`,f))
  T9<-data.table(
    quantity=c("Incoming influence (identity-swap null)",
               "Outgoing influence (timing-permutation null)",
               "Outgoing influence (circular-shift null)",
               "Outgoing influence (timing + marks null)"),
    real_median=round(c(med(ri),med(ro),med(ro),med(ro)),3),
    null_median=round(c(med(pool("inc_swap")),med(pool("out_perm")),
                        med(pool("out_circ")),med(pool("out_pm"))),3))
  print(T9); fwrite(T9,file.path(OUT_DIR,"table9.csv"))
  blk<-function(rv,fld,lbl){nl<-lapply(rows,`[[`,fld); k<-which(lengths(nl)>0)
  if(!length(k)) return(NULL)
  p<-mapply(emp_p,rv[k],nl[k]); nm<-vapply(nl[k],med,numeric(1))
  data.table(quantity=lbl,n_birds=length(k),
             real_med=round(med(rv[k]),4),null_med=round(med(nm),4),
             pct_above=round(100*mean(rv[k]>nm),1),
             pct_p_lt_05=round(100*mean(p<0.05),1),
             min_possible_p=round(1/(median(lengths(nl[k]))+1),4),
             fisher_p=signif(fisher(p),3))}
  T9s<-rbindlist(list(blk(ri,"inc_swap","Incoming (swap)"),
                      blk(ro,"out_perm","Outgoing (permutation)"),blk(ro,"out_circ","Outgoing (circular)"),
                      blk(ro,"out_pm","Outgoing (timing+marks)")))
  print(T9s); fwrite(T9s,file.path(OUT_DIR,"table9_stats.csv"))
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: bootstrap intervals (% of birds with CI excluding 0)

## --------------------------------------------------------------------
##  Section 7: block-restricted timing nulls
## --------------------------------------------------------------------
if (run("Sec7Blocks")) {
  message("\n### BLOCK-RESTRICTED TIMING NULLS ###")
  set.seed(SEED_BASE + 1110)
  
  ## K = 6 is omitted: the circular shift requires K*(2L+3) = 114 observed bins,
  ## which many flights do not have, so that row would be underpowered rather
  ## than informative. K = 1 and K = 3 carry the comparison.
  BN <- list(K = c(1, 3), N_PERM = 50)
  
  ## permute event times within contiguous blocks of the OBSERVED bins
  perm_blocks <- function(Yr, OBSr, K) {
    o <- which(OBSr == 1L); n <- length(o)
    if (n < 2*K) return(NULL)
    edges <- floor(seq(0, n, length.out = K + 1))
    out <- Yr
    for (k in seq_len(K)) {
      idx <- o[(edges[k]+1):edges[k+1]]
      if (length(idx) < 2) next
      out[idx] <- sample(Yr[idx])
    }
    out
  }
  
  ## circular shift within contiguous blocks; the wrap-around now falls at a
  ## block boundary rather than joining arrival to departure
  circ_blocks <- function(Yr, OBSr, K, L) {
    o <- which(OBSr == 1L); n <- length(o)
    if (n < K*(2*L + 3)) return(NULL)
    edges <- floor(seq(0, n, length.out = K + 1))
    out <- Yr
    for (k in seq_len(K)) {
      idx <- o[(edges[k]+1):edges[k+1]]; m <- length(idx)
      if (m < 2*L + 3) next
      s <- sample(seq(L+1, m-L-1), 1)
      out[idx] <- Yr[idx][((seq_len(m) - 1 + s) %% m) + 1]
    }
    out
  }
  
  ## outgoing null for source b, given a precomputed base design. Only column b
  ## of each target's design depends on b's events, so the rest is reused.
  out_null_row <- function(DZ, M, b, Yb, cfg, Sbase) {
    if (is.null(Yb)) return(NA_real_)
    tg <- setdiff(seq_len(DZ$m), b)
    So <- Sbase
    for (i in tg) So[[i]][, b] <- col_of(DZ, i, b, Yb, cfg)
    fit_flight(DZ, M, cfg, targets = tg, Sover = So)$out[b]
  }
  
  P <- prep_all()
  t0 <- Sys.time()
  for (idx in seq_along(P)) {
    fk <- names(P)[idx]
    ck <- file.path(OUT_DIR, "ckpt", paste0("bn_", gsub("[^A-Za-z0-9]","_",fk), ".rds"))
    if (file.exists(ck)) next
    G  <- P[[fk]]
    DZ <- build_W(G, CFG)
    f0 <- fit_flight(DZ, G, CFG)
    ## built once per flight, not once per null draw
    Sbase <- lapply(seq_len(DZ$m), function(i) design_i(DZ, G$Y, i, CFG))
    rows <- list()
    for (bi in seq_along(G$birds)) {
      rec <- list(fk = fk, bird = G$birds[bi], real_out = f0$out[bi])
      for (K in BN$K) {
        pn <- cn <- numeric(BN$N_PERM)
        for (r in seq_len(BN$N_PERM)) {
          pn[r] <- out_null_row(DZ, G, bi,
                                perm_blocks(G$Y[bi,], G$OBS[bi,], K), CFG, Sbase)
          cn[r] <- out_null_row(DZ, G, bi,
                                circ_blocks(G$Y[bi,], G$OBS[bi,], K, CFG$L), CFG, Sbase)
        }
        rec[[paste0("perm_K", K)]] <- pn[is.finite(pn)]
        rec[[paste0("circ_K", K)]] <- cn[is.finite(cn)]
      }
      rows[[length(rows)+1L]] <- rec
    }
    saveRDS(rows, ck)
    cat(sprintf("\r  %d/%d flights  %.0f min elapsed", idx, length(P),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    flush.console()
  }
  cat("\n")
  
  rows <- unlist(lapply(list.files(file.path(OUT_DIR,"ckpt"), pattern = "^bn_",
                                   full.names = TRUE), readRDS), recursive = FALSE)
  ro <- vapply(rows, `[[`, numeric(1), "real_out")
  
  blk <- function(fld, lbl, K) {
    nl <- lapply(rows, function(z) if (is.null(z[[fld]])) numeric(0) else z[[fld]])
    k <- which(lengths(nl) > 0); if (!length(k)) return(NULL)
    p  <- mapply(emp_p, ro[k], nl[k])
    nm <- vapply(nl[k], med, numeric(1))
    data.table(null = lbl, blocks = K, n_birds = length(k),
               real_med  = round(med(ro[k]), 4),
               null_med  = round(med(nm), 4),
               margin    = round(median(ro[k] - nm, na.rm = TRUE), 4),
               pct_above = round(100*mean(ro[k] > nm), 1),
               fisher_p  = signif(fisher(p), 3))
  }
  RES <- rbindlist(Filter(Negate(is.null), c(
    lapply(BN$K, function(K) blk(paste0("perm_K",K), "timing permutation", K)),
    lapply(BN$K, function(K) blk(paste0("circ_K",K), "circular shift", K)))))
  setorder(RES, null, blocks)
  print(RES); fwrite(RES, file.path(OUT_DIR, "blocknull_sensitivity.csv"))
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Section 7: null battery across event thresholds
## --------------------------------------------------------------------
if (run("Sec7Thresh")) {
  message("\n### DYADIC CLAIM ACROSS THRESHOLDS ###")
  set.seed(SEED_BASE + 1120)
  cfgD<-modifyList(CFG,list(N_PERM=25,N_SWAP=8))
  for(q in CFG$QUANTILES){
    qs<-gsub("\\.","",sprintf("%.2f",q)); Pq<-prep_all(q)
    bf<-list(); for(fk in names(Pq)) for(b in Pq[[fk]]$birds){
      k<-paste(flock_of(fk),b,sep="|"); bf[[k]]<-c(bf[[k]],fk)}
    for(fk in names(Pq)){
      ck<-file.path(OUT_DIR,"ckpt",sprintf("dy%s_%s.rds",qs,gsub("[^A-Za-z0-9]","_",fk)))
      if(file.exists(ck)) next
      G<-Pq[[fk]]; fl<-flock_of(fk); DZ<-build_W(G,cfgD); f0<-fit_flight(DZ,G,cfgD)
      rows<-list()
      for(bi in seq_along(G$birds)){
        b<-G$birds[bi]; dn<-setdiff(bf[[paste(fl,b,sep="|")]],fk); inc<-numeric(0)
        if(length(dn)){ if(length(dn)>cfgD$N_SWAP) dn<-sample(dn,cfgD$N_SWAP)
        for(d in dn){dr<-tryCatch(donor_reg(G,Pq[[d]]$reg,b,cfgD),error=function(e)NULL)
        if(is.null(dr)) next
        rs<-rbind(G$reg[bird!=b,.(bird,bin,x,y,disp,heading,event)],dr)
        Ms<-mats_from_reg(rs,G$birds,G$B)
        fs<-tryCatch(fit_flight(build_W(Ms,cfgD),Ms,cfgD),error=function(e)NULL)
        if(!is.null(fs)) inc<-c(inc,fs$inc[bi])}}
        pn<-cn<-numeric(cfgD$N_PERM)
        for(r in seq_len(cfgD$N_PERM)){pn[r]<-out_null(DZ,G,bi,"perm",cfgD)
        cn[r]<-out_null(DZ,G,bi,"circ",cfgD)}
        rows[[length(rows)+1L]]<-list(q=q,fk=fk,bird=b,real_inc=f0$inc[bi],
                                      real_out=f0$out[bi],inc_swap=inc,out_perm=pn[is.finite(pn)],
                                      out_circ=cn[is.finite(cn)])}
      saveRDS(rows,ck)}
    cat("q",q,"done\n")}
  DY<-list()
  for(q in CFG$QUANTILES){
    qs<-gsub("\\.","",sprintf("%.2f",q))
    fs<-list.files(file.path(OUT_DIR,"ckpt"),pattern=paste0("^dy",qs,"_"),full.names=TRUE)
    if(!length(fs)) next
    rows<-unlist(lapply(fs,readRDS),recursive=FALSE); if(!length(rows)) next
    ri<-vapply(rows,`[[`,numeric(1),"real_inc"); ro<-vapply(rows,`[[`,numeric(1),"real_out")
    b<-function(rv,fld,lbl){nl<-lapply(rows,`[[`,fld); k<-which(lengths(nl)>0)
    if(!length(k)) return(NULL)
    p<-mapply(emp_p,rv[k],nl[k]); nm<-vapply(nl[k],med,numeric(1))
    data.table(quantile=q,quantity=lbl,n_birds=length(k),
               real_med=round(med(rv[k]),4),null_med=round(med(nm),4),
               pct_above=round(100*mean(rv[k]>nm),1),fisher_p=signif(fisher(p),3))}
    DY[[paste0(q,"i")]]<-b(ri,"inc_swap","Incoming (swap)")
    DY[[paste0(q,"p")]]<-b(ro,"out_perm","Outgoing (permutation)")
    DY[[paste0(q,"c")]]<-b(ro,"out_circ","Outgoing (circular)")}
  DD<-rbindlist(Filter(Negate(is.null),DY)); print(DD)
  fwrite(DD,file.path(OUT_DIR,"dyadic_sensitivity.csv"))
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: residual diagnostics and Figure 3

## --------------------------------------------------------------------
##  Section 7: alternative event definitions
## --------------------------------------------------------------------
if (run("Sec7Events")) {
  message("\n### EVENT-DEFINITION SENSITIVITY ###")
  set.seed(SEED_BASE + 1130)
  need_pigeons()
  
  ## reference specification: displacement events at the median threshold.
  ## The directional-delay scores depend only on heading and are therefore
  ## the same for every event definition, which is what makes them a fixed
  ## external yardstick.
  REF <- DEL <- list()
  Pref <- prep_all(0.50, "displacement")
  for (fk in names(Pref)) {
    G <- Pref[[fk]]
    REF[[fk]] <- tryCatch(fit_flight(build_W(G, CFG), G, CFG, want_Ew = TRUE),
                          error = function(e) NULL)
    DEL[[fk]] <- tryCatch(delay_leader(list(Hd = G$Hd), CFG),
                          error = function(e) NULL)
  }
  
  TS <- list()
  for (ty in c("displacement", "acceleration", "turning")) for (q in CFG$QUANTILES) {
    Pq <- prep_all(q, ty); if (!length(Pq)) next
    fits <- list(); er <- c()
    for (fk in names(Pq)) {
      G <- Pq[[fk]]
      er <- c(er, sum(G$Y[G$OBS == 1L])/sum(G$OBS))
      f <- tryCatch(fit_flight(build_W(G, CFG), G, CFG, want_Ew = TRUE),
                    error = function(e) NULL)
      if (!is.null(f)) fits[[fk]] <- f
    }
    if (!length(fits)) next
    fv <- function(k) vapply(fits, `[[`, numeric(1), k)
    cm <- intersect(names(fits), names(REF))
    
    ## agreement of per-bird net influence with the reference specification
    rk <- vapply(cm, function(fk) {
      a <- fits[[fk]]$net; b <- REF[[fk]]$net
      if (length(a) != length(b) || sd(a) == 0 || sd(b) == 0) return(NA_real_)
      suppressWarnings(cor(a, b, method = "spearman"))
    }, numeric(1))
    
    ## agreement with the directional-delay benchmark, which does not depend
    ## on the event definition
    dl <- vapply(cm, function(fk) {
      d <- DEL[[fk]]; a <- fits[[fk]]$net
      if (is.null(d) || length(d) != length(a) || sd(d) == 0 || sd(a) == 0)
        return(NA_real_)
      suppressWarnings(cor(a, d, method = "spearman"))
    }, numeric(1))
    
    TS[[paste(ty, q)]] <- data.table(
      event_type = ty, quantile = q, flights = length(fits),
      event_rate = round(mean(er, na.rm = TRUE), 3),
      med_edges  = round(median(fv("n_edges")), 1),
      med_rho_R  = round(median(fv("rho_R")), 3),
      max_rho_R  = round(max(fv("rho_R")), 3),
      n_rho_ge1  = sum(fv("rho_R") >= 1),
      rank_agree_ref   = round(mean(rk, na.rm = TRUE), 3),
      rank_agree_delay = round(mean(dl, na.rm = TRUE), 3))
  }
  TT <- rbindlist(TS)
  print(TT); fwrite(TT, file.path(OUT_DIR, "threshold_sweep.csv"))
  
  message("  Complete")
}

## ---------------------------------------------------------------------
#  Generates Section 8.4: null battery repeated across event thresholds

## --------------------------------------------------------------------
##  Section 7: rho(R) against residual structure
## --------------------------------------------------------------------
if (run("Sec7RhoACF")) {
  message("\n### rho(R) AGAINST FLOCK-LEVEL RESIDUAL STRUCTURE ###")
  set.seed(SEED_BASE + 1140)
  
  P <- prep_all()
  
  ## helpers -------------------------------------------------------------
  lagcor <- function(a, b, k) {           # cor( a(t), b(t-k) ), NA-aware
    n <- length(a); if (k >= n) return(NA_real_)
    x <- a[(k+1):n]; y <- b[1:(n-k)]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 50) return(NA_real_)
    suppressWarnings(cor(x[ok], y[ok]))
  }
  cshift <- function(v, k) { n <- length(v); v[((seq_len(n) - 1 + k) %% n) + 1] }
  
  N_CTRL <- 20     # circular-shift replicates for the control
  
  RES <- list()
  for (fk in names(P)) {
    G  <- P[[fk]]
    f  <- tryCatch(fit_flight(build_W(G, CFG), G, CFG, want_fit = TRUE, want_Ew = TRUE),
                   error = function(e) NULL)
    if (is.null(f) || is.null(f$fitted)) next
    
    ## Bernoulli randomized quantile residuals on the common bin axis, so that
    ## a lag between two birds is a lag in real time
    ft <- f$fitted
    ft[, p_ := pmin(pmax(lambda, 1e-10), 1 - 1e-10)]
    ft[, r  := qnorm(pmin(pmax(runif(.N, ifelse(y == 0, 0, 1 - p_),
                                     ifelse(y == 0, 1 - p_, 1)), 1e-10), 1 - 1e-10))]
    Rm <- matrix(NA_real_, G$m, G$B)
    for (k in seq_len(nrow(ft))) Rm[ft$bird[k], ft$bin[k]] <- ft$r[k]
    ok <- which(rowSums(is.finite(Rm)) > 100)
    if (length(ok) < 2) next
    
    ## within-bird ACF at lag 5
    w5 <- vapply(ok, function(i) lagcor(Rm[i, ], Rm[i, ], 5), numeric(1))
    
    ## cross-bird, real
    c5 <- c(); c0 <- c()
    for (i in ok) for (j in ok) {
      if (i == j) next
      c5 <- c(c5, lagcor(Rm[i, ], Rm[j, ], 5))
      c0 <- c(c0, lagcor(Rm[i, ], Rm[j, ], 0))
    }
    
    ## cross-bird, circular-shift control: preserves each bird's own
    ## autocorrelation exactly, destroys alignment between birds
    k5 <- c(); k0 <- c()
    for (p in seq_len(N_CTRL)) {
      Rp <- Rm
      for (i in ok) Rp[i, ] <- cshift(Rm[i, ], sample(ncol(Rm), 1))
      for (i in ok) for (j in ok) {
        if (i == j) next
        k5 <- c(k5, lagcor(Rp[i, ], Rp[j, ], 5))
        k0 <- c(k0, lagcor(Rp[i, ], Rp[j, ], 0))
      }
    }
    
    RES[[fk]] <- data.table(
      frkey       = fk,
      flock       = sub("_.*", "", fk),
      m           = G$m,
      n_bins      = G$B,
      rho_R       = f$rho_R,
      rho_off     = f$rho_off,
      share_cross = 100 * f$shares[3],
      share_self  = 100 * f$shares[2],
      acf5_within = mean(w5, na.rm = TRUE),
      acf5_cross  = mean(c5, na.rm = TRUE),
      acf0_cross  = mean(c0, na.rm = TRUE),
      acf5_ctrl   = mean(k5, na.rm = TRUE),
      acf0_ctrl   = mean(k0, na.rm = TRUE))
    cat("."); flush.console()
  }
  cat("\n")
  D <- rbindlist(RES)
  D[, acf5_excess := acf5_cross - acf5_ctrl]   # shared component at lag 5
  D[, acf0_excess := acf0_cross - acf0_ctrl]   # shared component, contemporaneous
  fwrite(D, file.path(OUT_DIR, "rho_vs_acf_per_release.csv"))
  
  cat(sprintf("releases analysed: %d\n\n", nrow(D)))
  
  ## ---- descriptive -----------------------------------------------------
  cat("=== per-release distributions ===\n")
  print(D[, .(quantity = c("rho_R", "share_cross (%)", "acf5_within",
                           "acf5_excess", "acf0_excess"),
              median = round(c(median(rho_R), median(share_cross), median(acf5_within),
                               median(acf5_excess), median(acf0_excess)), 4),
              IQR_lo = round(c(quantile(rho_R,.25), quantile(share_cross,.25),
                               quantile(acf5_within,.25), quantile(acf5_excess,.25),
                               quantile(acf0_excess,.25)), 4),
              IQR_hi = round(c(quantile(rho_R,.75), quantile(share_cross,.75),
                               quantile(acf5_within,.75), quantile(acf5_excess,.75),
                               quantile(acf0_excess,.75)), 4))])
  
  ## ---- the association -------------------------------------------------
  ct <- function(x, y, lx, ly) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) return(NULL)
    s <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
    p <- suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
    data.table(x = lx, y = ly, n = sum(ok),
               spearman = round(unname(s$estimate), 3), sp_p = signif(s$p.value, 3),
               pearson  = round(unname(p$estimate), 3),
               pear_lo  = round(p$conf.int[1], 3), pear_hi = round(p$conf.int[2], 3))
  }
  A <- rbindlist(Filter(Negate(is.null), list(
    ct(D$rho_R,       D$acf5_excess, "rho_R",       "acf5_excess"),
    ct(D$rho_R,       D$acf0_excess, "rho_R",       "acf0_excess"),
    ct(D$rho_R,       D$acf5_within, "rho_R",       "acf5_within"),
    ct(D$share_cross, D$acf5_excess, "share_cross", "acf5_excess"),
    ct(D$share_cross, D$acf0_excess, "share_cross", "acf0_excess"),
    ct(D$share_cross, D$acf5_within, "share_cross", "acf5_within"),
    ct(D$share_self,  D$acf5_excess, "share_self",  "acf5_excess"))))
  cat("\n=== association across releases ===\n"); print(A)
  fwrite(A, file.path(OUT_DIR, "rho_vs_acf_correlations.csv"))
  
  ## ---- the split the decomposition sensitivity used --------------------
  D[, near_critical := rho_R > 0.95]
  cat("\n=== near-critical vs rest ===\n")
  print(D[, .(n = .N,
              rho_R       = round(median(rho_R), 3),
              share_cross = round(median(share_cross), 2),
              share_self  = round(median(share_self), 2),
              acf5_within = round(median(acf5_within), 4),
              acf5_excess = round(median(acf5_excess), 4),
              acf0_excess = round(median(acf0_excess), 4)), by = near_critical])
  for (v in c("acf5_excess", "acf0_excess", "acf5_within")) {
    w <- suppressWarnings(wilcox.test(D[[v]] ~ D$near_critical))
    cat(sprintf("  %s: Wilcoxon p = %.3g\n", v, w$p.value))
  }
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Table 10: group-activity decomposition
## --------------------------------------------------------------------
if (run("Table10")) {
  message("\n### FLOCK-LEVEL COMMON DRIVER ###")
  set.seed(SEED_BASE + 1150)
  
  H_WIN <- c(5, 10, 20)      # window for the flock term, in bins
  H_MAIN <- 10
  
  P <- prep_all()
  
  ## leave-one-out flock activity for target i: mean recent event rate of the
  ## OTHER birds. Predictable: uses bins b-1 .. b-H only.
  flock_z <- function(M, i, H) {
    m <- M$m; B <- M$B
    others <- setdiff(seq_len(m), i)
    if (!length(others)) return(rep(0, B))
    acc <- numeric(B)
    for (j in others) for (h in seq_len(H)) acc <- acc + lagi(M$Y[j, ], h, B)
    acc / (length(others) * H)
  }
  
  ## fit one flight with and without the flock term
  fit_with_flock <- function(G, cfg, H, use_flock = TRUE) {
    DZ <- build_W(G, cfg)
    m <- DZ$m
    alpha <- matrix(0, m, m); gam <- rep(NA_real_, m); mu <- rep(NA_real_, m)
    ll <- 0; np <- 0; sh <- matrix(NA_real_, m, 4)
    for (i in seq_len(m)) {
      at <- DZ$rows[[i]]; if (length(at) < 2) next
      y <- G$Y[i, at]; if (sum(y) < cfg$min_events) next
      S <- design_i(DZ, G$Y, i, cfg)[at, , drop = FALSE]
      Z <- if (use_flock) flock_z(G, i, H)[at] else NULL
      X <- if (use_flock) cbind(S, Z) else S
      keep <- which(colMeans(X) > cfg$min_colmass); if (!length(keep)) next
      A1 <- cbind(1, X[, keep, drop = FALSE]); p <- length(keep)
      nll <- function(th){lam<-pmax(as.numeric(A1%*%th),1e-12); -sum(y*log(lam)-lam)}
      gr  <- function(th){lam<-pmax(as.numeric(A1%*%th),1e-12)
      -as.numeric(crossprod(A1, y/lam-1))}
      o <- tryCatch(optim(c(max(mean(y),1e-6), rep(.1,p)), nll, gr, method="L-BFGS-B",
                          lower=c(cfg$eps_mu, rep(0,p)),
                          control=list(maxit=500, factr=1e10)), error=function(e) NULL)
      if (is.null(o)) next
      lam <- pmax(as.numeric(A1 %*% o$par), 1e-12)
      ll <- ll + sum(y*log(lam) - lam); np <- np + p + 1
      cf <- numeric(ncol(X)); cf[keep] <- o$par[-1]
      mu[i] <- o$par[1]
      alpha[i, ] <- cf[seq_len(m)]
      if (use_flock) gam[i] <- cf[m + 1]
      ## intensity decomposition: background / self / cross / flock
      con <- sweep(X[, keep, drop=FALSE], 2, o$par[-1], `*`)
      tot <- mean(o$par[1] + rowSums(con))
      if (is.finite(tot) && tot > 0) {
        selfcol  <- which(keep == i)
        flockcol <- if (use_flock) which(keep == m + 1) else integer(0)
        s_self  <- if (length(selfcol))  mean(con[, selfcol])  else 0
        s_flock <- if (length(flockcol)) mean(con[, flockcol]) else 0
        sh[i, ] <- c(o$par[1]/tot, s_self/tot,
                     (mean(rowSums(con)) - s_self - s_flock)/tot, s_flock/tot)
      }
    }
    off <- alpha; diag(off) <- 0
    list(alpha=alpha, off=off, gamma=gam, mu=mu, loglik=ll, npar=np,
         shares=colMeans(sh, na.rm=TRUE),
         rho_off=max(Mod(eigen(off, only.values=TRUE)$values)),
         n_edges=sum(off > cfg$tol_edge))
  }
  
  ## ---- main comparison at H = H_MAIN -----------------------------------
  CMP <- list()
  for (fk in names(P)) {
    G <- P[[fk]]
    f0 <- tryCatch(fit_with_flock(G, CFG, H_MAIN, use_flock=FALSE), error=function(e) NULL)
    f1 <- tryCatch(fit_with_flock(G, CFG, H_MAIN, use_flock=TRUE ), error=function(e) NULL)
    if (is.null(f0) || is.null(f1)) next
    CMP[[fk]] <- data.table(
      frkey=fk,
      cross_base = 100*f0$shares[3], cross_flock = 100*f1$shares[3],
      self_base  = 100*f0$shares[2], self_flock  = 100*f1$shares[2],
      mu_base    = 100*f0$shares[1], mu_flock    = 100*f1$shares[1],
      flock_share= 100*f1$shares[4],
      gamma_med  = median(f1$gamma, na.rm=TRUE),
      edges_base = f0$n_edges, edges_flock = f1$n_edges,
      ll_base=f0$loglik, ll_flock=f1$loglik,
      np_base=f0$npar,  np_flock=f1$npar)
    cat("."); flush.console()
  }
  cat("\n")
  D <- rbindlist(CMP)
  D[, cross_change := cross_flock - cross_base]
  fwrite(D, file.path(OUT_DIR,"flockterm_per_release.csv"))
  
  cat(sprintf("\nreleases compared: %d\n", nrow(D)))
  cat("\n=== INTENSITY DECOMPOSITION, POOLED ===\n")
  print(data.table(
    component = c("background","self","cross","flock term"),
    without   = round(c(mean(D$mu_base), mean(D$self_base), mean(D$cross_base), NA), 2),
    with      = round(c(mean(D$mu_flock), mean(D$self_flock), mean(D$cross_flock),
                        mean(D$flock_share)), 2)))
  
  cat("\n=== CHANGE IN CROSS-EXCITATION SHARE ===\n")
  cat(sprintf("  mean change   : %+.2f percentage points\n", mean(D$cross_change)))
  cat(sprintf("  median change : %+.2f\n", median(D$cross_change)))
  cat(sprintf("  releases where cross FELL: %.0f%%\n", 100*mean(D$cross_change < 0)))
  w <- suppressWarnings(wilcox.test(D$cross_flock, D$cross_base, paired=TRUE))
  cat(sprintf("  paired Wilcoxon p = %.3g\n", w$p.value))
  cat(sprintf("  relative reduction: %.1f%%\n",
              100*(mean(D$cross_base)-mean(D$cross_flock))/mean(D$cross_base)))
  
  dAIC <- (2*sum(D$np_base) - 2*sum(D$ll_base)) - (2*sum(D$np_flock) - 2*sum(D$ll_flock))
  cat(sprintf("\n  dAIC for adding the flock term: %.1f (positive favours the flock term)\n", dAIC))
  
  ## ---- does the reduction concentrate in the high-synchrony releases? ---
  if (file.exists(file.path(OUT_DIR,"rho_vs_acf_per_release.csv"))) {
    A <- fread(file.path(OUT_DIR,"rho_vs_acf_per_release.csv"))
    M2 <- merge(D, A[, .(frkey, acf0_excess, acf5_excess, rho_R)], by="frkey")
    cat("\n=== is the reduction largest where synchrony is largest? ===\n")
    cat(sprintf("  Spearman(acf0_excess, cross_change) = %.3f\n",
                suppressWarnings(cor(M2$acf0_excess, M2$cross_change, method="spearman",
                                     use="complete.obs"))))
    cat("  A strong NEGATIVE value means cross-excitation falls most in exactly the\n")
    cat("  releases with the most contemporaneous shared structure, which is what\n")
    cat("  the omitted-common-driver explanation predicts.\n")
    fwrite(M2, file.path(OUT_DIR,"flockterm_merged.csv"))
  }
  
  ## ---- window sensitivity ----------------------------------------------
  cat("\n=== WINDOW SENSITIVITY ===\n")
  SW <- list()
  for (H in H_WIN) {
    cb <- c(); cf <- c()
    for (fk in names(P)) {
      G <- P[[fk]]
      f0 <- tryCatch(fit_with_flock(G,CFG,H,FALSE), error=function(e) NULL)
      f1 <- tryCatch(fit_with_flock(G,CFG,H,TRUE ), error=function(e) NULL)
      if (is.null(f0)||is.null(f1)) next
      cb <- c(cb, 100*f0$shares[3]); cf <- c(cf, 100*f1$shares[3])
    }
    SW[[as.character(H)]] <- data.table(H_bins=H, H_sec=H*CFG$bin_sec,
                                        cross_base=round(mean(cb),2), cross_flock=round(mean(cf),2),
                                        reduction_pct=round(100*(mean(cb)-mean(cf))/mean(cb),1))
  }
  print(rbindlist(SW)); fwrite(rbindlist(SW), file.path(OUT_DIR,"flockterm_window.csv"))
  
  message("  Complete")
}

## --------------------------------------------------------------------
##  Section 7: null battery with the group-activity term
## --------------------------------------------------------------------
if (run("Sec7FlockNull")) {
  message("\n### NULL BATTERY UNDER THE FLOCK-TERM SPECIFICATION ###")
  set.seed(SEED_BASE + 1160)
  
  H_FLOCK <- 10                      # bins; 4 s, as in the flockterm block
  CKF <- file.path(OUT_DIR, "ckpt")
  
  ## leave-one-out flock activity, computed from a supplied event matrix so the
  ## nulls can pass in modified events
  flock_z_Y <- function(Y, i, H, B) {
    m <- nrow(Y); others <- setdiff(seq_len(m), i)
    if (!length(others)) return(rep(0, B))
    acc <- numeric(B)
    for (j in others) for (h in seq_len(H)) acc <- acc + lagi(Y[j, ], h, B)
    acc / (length(others) * H)
  }
  
  ## fit one flight, optionally with a supplied event matrix (for the nulls)
  fit_flock <- function(DZ, G, cfg, Ymat = NULL, targets = seq_len(DZ$m),
                        Sover = NULL) {
    Y <- if (is.null(Ymat)) G$Y else Ymat
    m <- DZ$m; alpha <- matrix(0, m, m); gam <- rep(NA_real_, m)
    for (i in targets) {
      at <- DZ$rows[[i]]; if (length(at) < 2) next
      y <- Y[i, at]; if (sum(y) < cfg$min_events) next
      S <- if (is.null(Sover)) design_i(DZ, Y, i, cfg) else Sover[[i]]
      Z <- flock_z_Y(Y, i, H_FLOCK, DZ$B)          # recomputed from THESE events
      X <- cbind(S, Z)[at, , drop = FALSE]
      keep <- which(colMeans(X) > cfg$min_colmass); if (!length(keep)) next
      A1 <- cbind(1, X[, keep, drop = FALSE]); p <- length(keep)
      nll <- function(th){lam<-pmax(as.numeric(A1%*%th),1e-12); -sum(y*log(lam)-lam)}
      gr  <- function(th){lam<-pmax(as.numeric(A1%*%th),1e-12)
      -as.numeric(crossprod(A1, y/lam-1))}
      o <- tryCatch(optim(c(max(mean(y),1e-6), rep(.1,p)), nll, gr, method="L-BFGS-B",
                          lower=c(cfg$eps_mu, rep(0,p)),
                          control=list(maxit=400, factr=1e10)), error=function(e) NULL)
      if (is.null(o)) next
      cf <- numeric(ncol(X)); cf[keep] <- o$par[-1]
      alpha[i, ] <- cf[seq_len(m)]; gam[i] <- cf[m + 1]
    }
    off <- alpha; diag(off) <- 0
    list(off=off, gamma=gam, out=colSums(off), inc=rowSums(off))
  }
  
  ## outgoing null: null source b's events, refit targets != b. The flock term is
  ## recomputed inside fit_flock from the modified event matrix.
  out_null_flock <- function(DZ, G, b, type, cfg) {
    Yb <- if (type == "perm") perm_ev(G$Y[b, ], G$OBS[b, ])
    else               shift_ev(G$Y[b, ], G$OBS[b, ], cfg$L)
    if (is.null(Yb)) return(NA_real_)
    Ymod <- G$Y; Ymod[b, ] <- Yb
    tg <- setdiff(seq_len(DZ$m), b)
    So <- lapply(seq_len(DZ$m), function(i) design_i(DZ, G$Y, i, cfg))
    for (i in tg) So[[i]][, b] <- col_of(DZ, i, b, Yb, cfg)
    fit_flock(DZ, G, cfg, Ymat = Ymod, targets = tg, Sover = So)$out[b]
  }
  
  P <- prep_all()
  bird_fl <- list()
  for (fk in names(P)) for (b in P[[fk]]$birds) {
    k <- paste(flock_of(fk), b, sep="|"); bird_fl[[k]] <- c(bird_fl[[k]], fk)
  }
  
  for (fk in names(P)) {
    ck <- file.path(CKF, paste0("t9f_", gsub("[^A-Za-z0-9]","_",fk), ".rds"))
    if (file.exists(ck)) next
    G <- P[[fk]]; fl <- flock_of(fk); DZ <- build_W(G, CFG)
    f0 <- fit_flock(DZ, G, CFG)
    rows <- list()
    for (bi in seq_along(G$birds)) {
      b <- G$birds[bi]
      ## incoming: identity swap
      dn <- setdiff(bird_fl[[paste(fl,b,sep="|")]], fk); inc <- numeric(0)
      if (length(dn)) {
        if (length(dn) > CFG$N_SWAP) dn <- sample(dn, CFG$N_SWAP)
        for (d in dn) {
          dr <- tryCatch(donor_reg(G, P[[d]]$reg, b, CFG), error=function(e) NULL)
          if (is.null(dr)) next
          rs <- rbind(G$reg[bird != b, .(bird,bin,x,y,disp,heading,event)], dr)
          Ms <- mats_from_reg(rs, G$birds, G$B)
          fs <- tryCatch(fit_flock(build_W(Ms,CFG), Ms, CFG), error=function(e) NULL)
          if (!is.null(fs)) inc <- c(inc, fs$inc[bi])
        }
      }
      ## outgoing: timing nulls
      pn <- cn <- numeric(CFG$N_PERM)
      for (r in seq_len(CFG$N_PERM)) {
        pn[r] <- out_null_flock(DZ, G, bi, "perm", CFG)
        cn[r] <- out_null_flock(DZ, G, bi, "circ", CFG)
      }
      rows[[length(rows)+1L]] <- list(fk=fk, bird=b,
                                      real_inc=f0$inc[bi], real_out=f0$out[bi], gamma=f0$gamma[bi],
                                      inc_swap=inc, out_perm=pn[is.finite(pn)], out_circ=cn[is.finite(cn)])
    }
    saveRDS(rows, ck); cat("."); flush.console()
  }
  cat("\n")
  
  rows <- unlist(lapply(list.files(CKF, pattern="^t9f_", full.names=TRUE), readRDS),
                 recursive = FALSE)
  ri <- vapply(rows,`[[`,numeric(1),"real_inc"); ro <- vapply(rows,`[[`,numeric(1),"real_out")
  pool <- function(f) unlist(lapply(rows,`[[`,f))
  
  T9F <- data.table(
    quantity = c("Incoming influence (identity-swap null)",
                 "Outgoing influence (timing-permutation null)",
                 "Outgoing influence (circular-shift null)"),
    real_median = round(c(med(ri), med(ro), med(ro)), 4),
    null_median = round(c(med(pool("inc_swap")), med(pool("out_perm")),
                          med(pool("out_circ"))), 4))
  print(T9F); fwrite(T9F, file.path(OUT_DIR,"table9_flockterm.csv"))
  
  blk <- function(rv, fld, lbl) {
    nl <- lapply(rows, `[[`, fld); k <- which(lengths(nl) > 0)
    if (!length(k)) return(NULL)
    p <- mapply(emp_p, rv[k], nl[k]); nm <- vapply(nl[k], med, numeric(1))
    data.table(quantity=lbl, n_birds=length(k),
               real_med=round(med(rv[k]),4), null_med=round(med(nm),4),
               diff_med=round(median(rv[k]-nm, na.rm=TRUE),4),
               pct_above=round(100*mean(rv[k] > nm),1),
               pct_p_lt_05=round(100*mean(p < 0.05),1),
               min_possible_p=round(1/(median(lengths(nl[k]))+1),4),
               fisher_p=signif(fisher(p),3))
  }
  T9FS <- rbindlist(list(
    blk(ri,"inc_swap","Incoming (swap)"),
    blk(ro,"out_perm","Outgoing (permutation)"),
    blk(ro,"out_circ","Outgoing (circular)")))
  print(T9FS); fwrite(T9FS, file.path(OUT_DIR,"table9_flockterm_stats.csv"))
  
  cat(sprintf("\nflock-term loading gamma: median %.4f, nonzero in %.0f%% of birds\n",
              med(pool("gamma")), 100*mean(pool("gamma") > CFG$tol_edge, na.rm=TRUE)))
  
  ## side-by-side with the original specification, if present
  if (file.exists(file.path(OUT_DIR,"table9_stats.csv"))) {
    old <- fread(file.path(OUT_DIR,"table9_stats.csv"))
    cat("\n=== BASELINE vs FLOCK-TERM SPECIFICATION ===\n")
    cat("baseline columns:", paste(names(old), collapse=", "), "\n")
    pick <- function(dt, want) {          # first matching name, else NA
      hit <- want[want %in% names(dt)]
      if (length(hit)) dt[[hit[1]]] else rep(NA_real_, nrow(dt))
    }
    base <- data.table(
      quantity    = pick(old, c("quantity")),
      base_real   = pick(old, c("real_med","real_median")),
      base_null   = pick(old, c("null_med","null_median")),
      base_pct    = pick(old, c("pct_above","pct_above_null")),
      base_fisher = pick(old, c("fisher_p")))
    cmp <- merge(base,
                 T9FS[, .(quantity, flock_real = real_med, flock_null = null_med,
                          flock_pct = pct_above, flock_fisher = fisher_p)],
                 by = "quantity", all = TRUE)
    print(cmp); fwrite(cmp, file.path(OUT_DIR,"table9_comparison.csv"))
  }
  message("  Complete")
}

## --------------------------------------------------------------------
##  Figure 1: representative flight
## --------------------------------------------------------------------
if (run("Figure1")) {
  message("\n### FIGURE 1: REPRESENTATIVE FLIGHT ###")
  set.seed(SEED_BASE + 1170)
  
  FIG <- list(flight  = NULL,   # NULL selects the flight with median rho(R)
              R_BOOT  = 500,
              lwd_max = 7,
              node_r  = 0.155,  # node radius in plot units
              curv    = 0.20)   # arc curvature
  
  P <- prep_all()
  
  ## representative flight: rho(R) closest to the median across releases
  if (is.null(FIG$flight)) {
    rr <- vapply(names(P), function(fk) {
      f <- tryCatch(fit_flight(build_W(P[[fk]], CFG), P[[fk]], CFG, want_Ew = TRUE),
                    error = function(e) NULL)
      if (is.null(f)) NA_real_ else f$rho_R
    }, numeric(1))
    rr <- rr[is.finite(rr)]
    FIG$flight <- names(rr)[which.min(abs(rr - median(rr)))]
    cat(sprintf("  representative flight: %s (rho(R) = %.3f, median %.3f)\n",
                FIG$flight, rr[FIG$flight], median(rr)))
  }
  
  G  <- P[[FIG$flight]]
  DZ <- build_W(G, CFG)
  f0 <- fit_flight(DZ, G, CFG, want_Ew = TRUE)
  
  ## bootstrap over bins within each target row, to select edges
  Sall <- lapply(seq_len(DZ$m), function(i) design_i(DZ, G$Y, i, CFG))
  LO <- matrix(NA_real_, DZ$m, DZ$m)
  BT <- array(NA_real_, c(FIG$R_BOOT, DZ$m, DZ$m))
  for (r in seq_len(FIG$R_BOOT)) {
    for (i in seq_len(DZ$m)) {
      at <- DZ$rows[[i]]; if (length(at) < 2) next
      idx <- sample(at, length(at), replace = TRUE)
      fr <- fit_row(Sall[[i]][idx, , drop = FALSE], G$Y[i, idx], CFG)
      if (!is.null(fr)) BT[r, i, ] <- fr$alpha
    }
  }
  for (i in seq_len(DZ$m)) for (j in seq_len(DZ$m))
    LO[i, j] <- quantile(BT[, i, j], 0.025, na.rm = TRUE, names = FALSE)
  
  off <- f0$alpha; diag(off) <- 0
  sel <- (off > CFG$tol_edge) & (LO > 0); diag(sel) <- FALSE
  net <- f0$net
  
  EDG <- data.table(
    source = G$birds[col(off)[sel]], target = G$birds[row(off)[sel]],
    alpha  = round(off[sel], 4), lower = round(LO[sel], 4))
  setorder(EDG, -alpha)
  cat(sprintf("  %d of %d directed edges selected\n", sum(sel), DZ$m*(DZ$m-1)))
  print(EDG)
  fwrite(EDG, file.path(OUT_DIR, "figure1_edges.csv"))
  
  SCORE <- data.table(bird = G$birds, net = round(net, 4),
                      out = round(f0$out, 4), inc = round(f0$inc, 4),
                      self = round(diag(f0$alpha), 4))
  SCORE[, role := fifelse(bird == G$birds[which.max(net)], "Top leader",
                          fifelse(net > 0, "Leader/initiator", "Follower"))]
  print(SCORE); fwrite(SCORE, file.path(OUT_DIR, "figure1_scores.csv"))
  
  ## ---- drawing helpers ------------------------------------------------
  ramp <- colorRamp(c("#3B6DB5", "#F2F2F2", "#C0392B"))
  nlim <- max(abs(net), na.rm = TRUE); nlim <- if (nlim > 0) nlim else 1
  nodecol <- function(v) {
    z <- (v + nlim)/(2*nlim); z <- pmin(pmax(z, 0), 1)
    rgb(ramp(z), maxColorValue = 255)
  }
  circ <- function(x, y, r, ...) {
    th <- seq(0, 2*pi, length.out = 80)
    polygon(x + r*cos(th), y + r*sin(th), ...)
  }
  ## quadratic Bezier from (x0,y0) to (x1,y1), trimmed at both node boundaries
  arc <- function(x0, y0, x1, y1, w, col, curv, r) {
    dx <- x1 - x0; dy <- y1 - y0; L <- sqrt(dx^2 + dy^2)
    if (L < 1e-8) return(invisible())
    ux <- dx/L; uy <- dy/L
    x0 <- x0 + ux*r; y0 <- y0 + uy*r
    x1 <- x1 - ux*r; y1 <- y1 - uy*r
    mx <- (x0 + x1)/2; my <- (y0 + y1)/2
    cx <- mx - uy*curv*L; cy <- my + ux*curv*L
    t  <- seq(0, 1, length.out = 60)
    bx <- (1-t)^2*x0 + 2*(1-t)*t*cx + t^2*x1
    by <- (1-t)^2*y0 + 2*(1-t)*t*cy + t^2*y1
    lines(bx, by, lwd = w, col = col, lend = 1)
    n <- length(bx)
    arrows(bx[n-3], by[n-3], bx[n], by[n], length = 0.10,
           angle = 22, lwd = w, col = col)
  }
  
  ## ---- draw -----------------------------------------------------------
  png(file.path(OUT_DIR, "pigeon_representative_flight.png"),
      width = 2600, height = 1150, res = 210)
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 0.30, 1.05))
  
  ## ---- panel (a): influence network -----------------------------------
  par(mar = c(4.2, 4.2, 4.0, 0.5), xpd = FALSE)
  m   <- DZ$m
  ord <- order(net, decreasing = TRUE)
  ang <- pi/2 - seq(0, 2*pi, length.out = m + 1)[-(m+1)]   # top, then clockwise
  px  <- cos(ang); py <- sin(ang)
  nx  <- numeric(m); ny <- numeric(m)
  nx[ord] <- px; ny[ord] <- py
  
  plot(NA, xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4), asp = 1,
       axes = FALSE, xlab = "", ylab = "")
  par(xpd = NA)                      # labels and title sit outside the region
  title(main = "(a) Estimated influence network", adj = 0, cex.main = 1.15)
  mtext(expression("arrow " * j %->% i * " weighted by " * hat(alpha)[ij] *
                     "; node shaded by net leadership"),
        side = 3, line = 0.3, adj = 0, cex = 0.72)
  
  if (nrow(EDG)) {
    amax <- max(EDG$alpha)
    for (k in rev(seq_len(nrow(EDG)))) {   # weakest first, strong on top
      a <- which(G$birds == EDG$source[k]); b <- which(G$birds == EDG$target[k])
      w <- 1 + FIG$lwd_max*EDG$alpha[k]/amax
      shade <- grey(0.72 - 0.42*EDG$alpha[k]/amax)
      arc(nx[a], ny[a], nx[b], ny[b], w, shade, FIG$curv, FIG$node_r)
    }
  }
  for (i in seq_len(m)) {
    circ(nx[i], ny[i], FIG$node_r, col = nodecol(net[i]),
         border = grey(0.25), lwd = 1.4)
    text(nx[i], ny[i], G$birds[i], cex = 0.82, font = 2)
    text(nx[i]*1.32, ny[i]*1.32, sprintf("%+.3f", net[i]),
         cex = 0.70, col = grey(0.30))
  }
  
  ## ---- colour bar ------------------------------------------------------
  par(mar = c(6, 0.5, 6, 3.0), xpd = NA)
  plot(NA, xlim = c(0, 1), ylim = c(-nlim, nlim), axes = FALSE,
       xlab = "", ylab = "")
  nyb <- 200
  yy  <- seq(-nlim, nlim, length.out = nyb + 1)
  for (k in seq_len(nyb))
    rect(0.05, yy[k], 0.55, yy[k+1], col = nodecol((yy[k] + yy[k+1])/2),
         border = NA)
  rect(0.05, -nlim, 0.55, nlim, border = grey(0.35))
  tk <- pretty(c(-nlim, nlim), 5); tk <- tk[tk >= -nlim & tk <= nlim]
  segments(0.55, tk, 0.62, tk, col = grey(0.35))
  text(0.68, tk, sprintf("%+.2f", tk), adj = 0, cex = 0.75)
  text(0.05, nlim*1.16, "Net\nleadership", adj = 0, cex = 0.80, font = 2)
  
  ## ---- panel (b): flight tracks by recovered role -----------------------
  ## xpd = FALSE so that grid() is clipped to this panel
  par(mar = c(4.2, 4.6, 4.0, 1.0), xpd = FALSE)
  D <- PIG[frkey == FIG$flight]
  rolecol <- c("Follower" = "#3B6DB5", "Leader/initiator" = "#E8913A",
               "Top leader" = "#C0392B")
  plot(NA, xlim = range(D$lon, na.rm = TRUE), ylim = range(D$lat, na.rm = TRUE),
       xlab = "Longitude", ylab = "Latitude")
  grid(col = grey(0.90), lty = 1)
  for (rl in c("Follower", "Leader/initiator", "Top leader")) {
    for (b in SCORE$bird[SCORE$role == rl]) {
      d <- D[bird == b]
      lines(d$lon, d$lat, col = rolecol[rl], lwd = 1.8)
      if (nrow(d)) points(d$lon[nrow(d)], d$lat[nrow(d)], pch = 16,
                          col = rolecol[rl], cex = 0.9)
    }
  }
  legend("bottomleft", legend = names(rolecol), col = rolecol, lwd = 2.2,
         pch = 16, bty = "n", cex = 0.80, inset = c(0.02, 0.02))
  par(xpd = NA)                      # title only
  title(main = "(b) Flight tracks by recovered role", adj = 0, cex.main = 1.15)
  mtext(sprintf("%s; top leader = %s", FIG$flight, G$birds[which.max(net)]),
        side = 3, line = 0.3, adj = 0, cex = 0.72)
  dev.off()
  
  cat("  wrote pigeon_representative_flight.png, figure1_edges.csv and figure1_scores.csv\n")
  cat("\n  Arrow widths are proportional to the fitted alpha. Only edges whose\n")
  cat("  bootstrap lower bound excludes zero are drawn, which is the selection\n")
  cat("  rule evaluated in Section 6 rather than an arbitrary threshold. Net\n")
  cat("  leadership is printed beside each node and encoded in its fill.\n")
  message("  Complete")
}

message("\nDone. Outputs written to ", normalizePath(OUT_DIR))
sessionInfo()
