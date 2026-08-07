#!/usr/bin/env Rscript
## =====================================================================
##  hawkes_trajectories_analysis.R
##
##  Master script. Regenerates every numbered table and every reported
##  quantity in "Hawkes Trajectories: A Framework for Leadership Inference
##  in Coordinated Movement".
##
##  USAGE
##    Rscript hawkes_trajectories_analysis.R              # everything
##    Rscript hawkes_trajectories_analysis.R sim          # simulations only
##    Rscript hawkes_trajectories_analysis.R table9 diag  # named blocks
##
##  BLOCKS: sim table9 realfits boot thresh dyadic diag lag5 origin marks
##
##  Pigeon data (Watts et al. 2016, Dryad doi:10.5061/dryad.508j2) in ./Data.
##  Dependencies: data.table, sf.
## =====================================================================

suppressPackageStartupMessages({ library(data.table); library(sf) })
set.seed(20260806)

ARGS   <- commandArgs(trailingOnly = TRUE)
BLOCKS <- if (length(ARGS)) ARGS else
  c("sim","table9","realfits","boot","thresh","dyadic","diag","lag5","origin","marks")
run <- function(b) b %in% BLOCKS
DATA_DIR <- "Data"; OUT_DIR <- "."
dir.create(file.path(OUT_DIR,"ckpt"), showWarnings=FALSE, recursive=TRUE)

CFG <- list(
  pattern="\\.csv$", utm_epsg=32630,
  speed_kmh_fly=15, min_fly_sec=30,
  bin_sec=0.4, L=8,
  beta=3.0, kappa=6, sigma_s=15, sigma_d=15, d_max=60,
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
        w<-wh_of(angd(tg[[i]]$th,hj),cfg$kappa)*wd_of(d2,cfg$sigma_s,cfg$d_max)
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
## ---------------------------------------------------------------------
if (run("sim")) {
message("\n### SIMULATION STUDY ###")
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

# Generates Table 7 (Study 5, bin width and working likelihood)
S5<-list(); At<-make_A("single",CFG$M)
for(bs in c(0.2,0.4,0.8,1.6)) for(fam in c("poisson","bernoulli")){
  c2<-CFG; c2$bin_sec<-bs; c2$B<-round(CFG$B*CFG$bin_sec/bs)
  er<-c();rc<-c();evr<-c();k<-0
  for(r in seq_len(60)){
    mk<-simulate_marks(c2); W<-build_W_sim(mk,eta0,c2)
    Y<-simulate_Y(W,At,c2$mu,c2,fam); if(min(rowSums(Y))<10) next
    f<-fit_sim(W,Y,c2); if(is.null(f)) next
    k<-k+1; m<-metrics(f$A,At,c2$tol_edge)
    er<-c(er,m$edge_rmse);rc<-c(rc,m$rank_cor)
    evr<-c(evr,mean(Y[,(c2$L+1):c2$B]>0))}
  if(!k) next
  S5[[paste(bs,fam)]]<-data.table(bin_sec=bs,family=fam,replicates=k,
    event_rate=round(mean(evr),3),edge_RMSE=round(mean(er),4),
    edge_RMSE_mcse=round(mcse(er),5),rank_cor=round(mean(rc,na.rm=TRUE),3),
    rank_max=round(rank_ceiling(net_lead(At)),3))}
T7<-rbindlist(S5); print(T7); fwrite(T7,file.path(OUT_DIR,"table_sim_binwidth.csv"))
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: rho(R), intensity shares, dAIC for free diagonal
## ---------------------------------------------------------------------
if (run("realfits")) {
message("\n### PIGEON REAL FITS ###")
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
}

## ---------------------------------------------------------------------
#  Generates Table 9 (null battery) and Section 8.3 p-values
## ---------------------------------------------------------------------
if (run("table9")) {
message("\n### TABLE 9: NULL BATTERY ###")
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
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: bootstrap intervals (% of birds with CI excluding 0)
## ---------------------------------------------------------------------
if (run("boot")) {
message("\n### CONSTRAINED PARAMETRIC BOOTSTRAP ###")
P<-prep_all()
sim_Y <- function(DZ,M,mu,alpha,cfg){
  m<-DZ$m;B<-DZ$B; Ys<-matrix(0,m,B); Ys[,1:cfg$L]<-M$Y[,1:cfg$L]
  for(b in (cfg$L+1):B) for(i in seq_len(m)){
    if(M$OBS[i,b]!=1L){Ys[i,b]<-0;next}
    lam<-mu[i]
    for(j in seq_len(m)) for(l in seq_len(cfg$L)){
      bl<-b-l; if(bl<1) next
      lam<-lam+alpha[i,j]*DZ$W[[i]][[j]][[l]][b]*Ys[j,bl]}
    Ys[i,b]<-rbinom(1,1,min(max(lam,0),1))}
  Ys}
BR<-list()
for(fk in names(P)){
  ck<-file.path(OUT_DIR,"ckpt",paste0("bt_",gsub("[^A-Za-z0-9]","_",fk),".rds"))
  if(file.exists(ck)){BR[[fk]]<-readRDS(ck); next}
  G<-P[[fk]]; DZ<-build_W(G,CFG); f0<-fit_flight(DZ,G,CFG)
  mu<-f0$mu; mu[!is.finite(mu)]<-CFG$eps_mu
  bo<-vector("list",CFG$R_BOOT)
  for(r in seq_len(CFG$R_BOOT)){
    Ys<-sim_Y(DZ,G,mu,f0$alpha,CFG); Mb<-G; Mb$Y<-Ys
    fb<-tryCatch(fit_flight(DZ,Mb,CFG),error=function(e)NULL)
    if(!is.null(fb)) bo[[r]]<-list(out=fb$out,inc=fb$inc,net=fb$net)}
  bo<-Filter(Negate(is.null),bo)
  res<-list(birds=G$birds,hat=f0,boot=bo); saveRDS(res,ck); BR[[fk]]<-res; cat(".")}
cat("\n")
BD<-rbindlist(lapply(names(BR),function(fk){z<-BR[[fk]]; if(!length(z$boot)) return(NULL)
  rbindlist(lapply(seq_along(z$birds),function(k){
    q<-function(f) quantile(vapply(z$boot,function(b) b[[f]][k],numeric(1)),
                            c(.025,.975),names=FALSE,na.rm=TRUE)
    qo<-q("out");qi<-q("inc");qn<-q("net")
    data.table(fk=fk,bird=z$birds[k],out=z$hat$out[k],out_lo=qo[1],out_hi=qo[2],
      inc=z$hat$inc[k],inc_lo=qi[1],inc_hi=qi[2],
      net=z$hat$net[k],net_lo=qn[1],net_hi=qn[2])}))}))
cat(sprintf("birds with influence CI excluding zero: %.1f%%\n",
            100*mean(BD$out_lo>0|BD$inc_lo>0,na.rm=TRUE)))
fwrite(BD,file.path(OUT_DIR,"bootstrap_intervals.csv"))
}

## ---------------------------------------------------------------------
#  Generates Section 8.4: event-definition sensitivity (thresholds and types)
## ---------------------------------------------------------------------
if (run("thresh")) {
message("\n### EVENT-DEFINITION SENSITIVITY ###")
need_pigeons()
REF<-DEL<-list()
Pref<-prep_all(0.50,"displacement")
for(fk in names(Pref)){G<-Pref[[fk]]
  REF[[fk]]<-tryCatch(fit_flight(build_W(G,CFG),G,CFG),error=function(e)NULL)
  DEL[[fk]]<-tryCatch(delay_leader(list(Hd=G$Hd),CFG),error=function(e)NULL)}
TS<-list()
for(ty in c("displacement","acceleration","turning")) for(q in CFG$QUANTILES){
  Pq<-prep_all(q,ty); if(!length(Pq)) next
  fits<-list(); er<-c()
  for(fk in names(Pq)){G<-Pq[[fk]]
    er<-c(er,sum(G$Y[G$OBS==1L])/sum(G$OBS))
    f<-tryCatch(fit_flight(build_W(G,CFG),G,CFG),error=function(e)NULL)
    if(!is.null(f)) fits[[fk]]<-f}
  if(!length(fits)) next
  fv<-function(k) vapply(fits,`[[`,numeric(1),k)
  cm<-intersect(names(fits),names(REF))
  rk<-vapply(cm,function(fk){a<-fits[[fk]]$net;b<-REF[[fk]]$net
    if(length(a)!=length(b)||sd(a)==0||sd(b)==0) return(NA_real_)
    suppressWarnings(cor(a,b,method="spearman"))},numeric(1))
  dl<-vapply(cm,function(fk){d<-DEL[[fk]];a<-fits[[fk]]$net
    if(is.null(d)||length(d)!=length(a)||sd(d)==0||sd(a)==0) return(NA_real_)
    suppressWarnings(cor(a,d,method="spearman"))},numeric(1))
  TS[[paste(ty,q)]]<-data.table(event_type=ty,quantile=q,flights=length(fits),
    event_rate=round(mean(er,na.rm=TRUE),3),med_edges=round(median(fv("n_edges")),1),
    med_rho_R=round(median(fv("rho_R")),3),max_rho_R=round(max(fv("rho_R")),3),
    rank_agree_ref=round(mean(rk,na.rm=TRUE),3),
    rank_agree_delay=round(mean(dl,na.rm=TRUE),3))}
TT<-rbindlist(TS); print(TT); fwrite(TT,file.path(OUT_DIR,"threshold_sweep.csv"))
}

## ---------------------------------------------------------------------
#  Generates Section 8.4: null battery repeated across event thresholds
## ---------------------------------------------------------------------
if (run("dyadic")) {
message("\n### DYADIC CLAIM ACROSS THRESHOLDS ###")
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
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: residual diagnostics and Figure 3
## ---------------------------------------------------------------------
if (run("diag")) {
message("\n### RESIDUAL DIAGNOSTICS ###")
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
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: cross-bird check on the 2 s residual oscillation
## ---------------------------------------------------------------------
if (run("lag5")) {
message("\n### 2 s OSCILLATION: WITHIN- AND CROSS-BIRD ###")
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
}

## ---------------------------------------------------------------------
#  Generates Section 8.2: origin tests for the 2 s oscillation
## ---------------------------------------------------------------------
if (run("origin")) {
message("\n### 2 s OSCILLATION: ORIGIN ###")
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
}

## ---------------------------------------------------------------------
#  Generates Section 6 note: predictable vs contemporaneous target marks
## ---------------------------------------------------------------------
if (run("marks")) {
message("\n### PREDICTABLE-MARKS SENSITIVITY ###")
P<-prep_all(); LAG<-CON<-list()
for(fk in names(P)){G<-P[[fk]]
  LAG[[fk]]<-tryCatch(fit_flight(build_W(G,CFG,TRUE ),G,CFG),error=function(e)NULL)
  CON[[fk]]<-tryCatch(fit_flight(build_W(G,CFG,FALSE),G,CFG),error=function(e)NULL)}
LAG<-Filter(Negate(is.null),LAG); CON<-Filter(Negate(is.null),CON)
fv<-function(L,k) vapply(L,`[[`,numeric(1),k)
sh<-function(L) rowMeans(vapply(L,`[[`,numeric(3),"shares"),na.rm=TRUE)
ov<-function(L) unlist(lapply(L,function(z) z$off[z$off>CFG$tol_edge]))
sL<-sh(LAG); sC<-sh(CON)
MK<-data.table(quantity=c("median off-diagonal alpha","median rho(R)","median alpha_ii",
  "background %","self %","cross %","total loglik"),
  lagged=c(round(median(ov(LAG)),4),round(median(fv(LAG,"rho_R")),4),
    round(median(unlist(lapply(LAG,`[[`,"self")),na.rm=TRUE),4),
    round(100*sL[1],1),round(100*sL[2],1),round(100*sL[3],1),round(sum(fv(LAG,"loglik")),1)),
  contemporaneous=c(round(median(ov(CON)),4),round(median(fv(CON,"rho_R")),4),
    round(median(unlist(lapply(CON,`[[`,"self")),na.rm=TRUE),4),
    round(100*sC[1],1),round(100*sC[2],1),round(100*sC[3],1),round(sum(fv(CON,"loglik")),1)))
print(MK); fwrite(MK,file.path(OUT_DIR,"marks_sensitivity.csv"))
cat(sprintf("dAIC favouring lagged: %.1f\n",
            2*(sum(fv(LAG,"loglik"))-sum(fv(CON,"loglik")))))
}

message("\nDone. Outputs written to ", normalizePath(OUT_DIR))
sessionInfo()

