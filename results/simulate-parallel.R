rm(list=ls())

args <- commandArgs()
print(args)

library(data.table)
library(magrittr)
library(R.utils)
library(stringr)
library(flexsurv)

setwd("~/repos/XSRecency/")
source("./results/sim-helpers.R")
source("./R/phi-functions.R")

# Get command-line arguments
a <- commandArgs(trailingOnly=TRUE, asValues=TRUE,
                    defaults=list(
                      seed=100,
                      n_sims=2,
                      n=1000,
                      p=0.2,
                      inc=0.03,
                      window=101,
                      shadow=194,
                      itype="constant",
                      rho=0.0,
                      tau=10,
                      bigT=1,
                      phi_frr=NULL,
                      phi_tfrr=NULL,
                      phi_norm_mu=NULL,
                      phi_norm_sd=NULL,
                      phi_norm_div=NULL,
                      phi_pnorm_mu=NULL,
                      phi_pnorm_sd=NULL,
                      phi_pnorm_div=NULL,
                      out_dir=".",
                      ext_FRR=FALSE,
                      duong_scale=NULL,
                      max_FRR=NULL,
                      last_point=FALSE,
                      pt=TRUE,
                      t_min=0,
                      t_max=1,
                      q=1,
                      gamma=0, # variance for the Gaussian noise to add to prior test time
                      eta=0, # the probability of incorrectly reporting negative test
                      nu=0, # the probability of failing to report prior test result
                      xi=0, # the probability of failing to report prior positive test results
                      mech2=FALSE
                    ))

# Capture date in the out directory
date <- format(Sys.time(), "%d-%m-%y-%H")
out_dir <- paste0(a$out_dir, "/", date, "/")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

a$out_dir <- NULL

if(!is.null(a$rho)) a$rho <- as.numeric(a$rho)
if(!is.null(a$phi_frr)) a$phi_frr <- as.numeric(a$phi_frr)
if(!is.null(a$phi_tfrr)) a$phi_tfrr <- as.numeric(a$phi_tfrr)
if(!is.null(a$phi_norm_mu)) a$phi_norm_mu <- as.numeric(a$phi_norm_mu)
if(!is.null(a$phi_norm_sd)) a$phi_norm_sd <- as.numeric(a$phi_norm_sd)
if(!is.null(a$phi_norm_div)) a$phi_norm_div <- as.numeric(a$phi_norm_div)
if(!is.null(a$phi_pnorm_mu)) a$phi_pnorm_mu <- as.numeric(a$phi_pnorm_mu)
if(!is.null(a$phi_pnorm_sd)) a$phi_pnorm_sd <- as.numeric(a$phi_pnorm_sd)
if(!is.null(a$phi_pnorm_div)) a$phi_pnorm_div <- as.numeric(a$phi_pnorm_div)
if(!is.null(a$max_FRR)) a$max_FRR <- as.numeric(a$max_FRR)
if(!is.null(a$t_min)) a$t_min <- as.numeric(a$t_min)
if(!is.null(a$t_max)) a$t_max <- as.numeric(a$t_max)
if(!is.null(a$q)) a$q <- as.numeric(a$q)
if(!is.null(a$gamma)) a$gamma <- as.numeric(a$gamma)
if(!is.null(a$eta)) a$eta <- as.numeric(a$eta)
if(!is.null(a$nu)) a$nu <- as.numeric(a$nu)
if(!is.null(a$xi)) a$xi <- as.numeric(a$xi)

# Logic checks for arguments
if(!is.null(a$phi_frr) & !is.null(a$phi_tfrr)){
  stop("Can't provide both frr and time for frr.")
}
if(!is.null(a$phi_norm_mu)){
  if(is.null(a$phi_norm_sd) | is.null(a$phi_norm_div)){
    stop("Need stdev and divided by params for normal.")
  }
}
if(is.null(a$rho) & a$itype != "constant"){
  stop("Need a rho param if not constant incidence.")
}
if(is.null(a$rho)) a$rho <- NA

# Get the gamma parameters and baseline phi function
params <- get.gamma.params(window=a$window/365.25, shadow=a$shadow/365.25)

# Set up each type of phi function, will be overwritten
phi.none <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
phi.const <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
phi.norm <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
phit.pnorm <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

phi.func <- phi.none

# Get the phi function with constant FRR either past a certain time
# or fixed after it hits some value.
if(!is.null(a$phi_tfrr) | !is.null(a$phi_frr)){
  if(!is.null(a$phi_tfrr)){
    ttime <- a$phi_tfrr
    tval <- phi.none(ttime)
  }
  if(!is.null(a$phi_frr)){
    tval <- a$phi_frr
    ttime <- uniroot(function(t) phi.none(t) - tval, interval=c(0, a$tau))$root
  }
  phi.const <- function(t, ...) phi.none(t)*(t <= ttime) + tval*(t > ttime)
  phi.func <- phi.const
}
if(!is.null(a$phi_norm_mu)){
  phi.norm <- function(t) phi.const(t) + dnorm(t-a$phi_norm_mu, mean=0, sd=a$phi_norm_sd) / a$phi_norm_div
  phi.func <- phi.norm
}
if(!is.null(a$phi_pnorm_mu)){
  phi.pnorm <- function(t) phi.const(t) + pnorm(t-a$phi_pnorm_mu, mean=0, sd=a$phi_pnorm_sd) / a$phi_pnorm_div
  phi.func <- phi.pnorm
}

if(a$itype == "constant"){
  inc.function <- incidence.con
  infection.function <- infections.con
} else if(a$itype == "linear"){
  inc.function <- incidence.lin
  infection.function <- infections.lin
} else if(a$itype == "exponential"){
  inc.function <- incidence.exp
  infection.function <- infections.exp
} else if(a$itype == "piecewise"){
  # TODO: Figure out how to generate infection times from this.
  inc.function <- NULL
  infection.funtion <- NULL
} else {
  stop("Unknown incidence function.")
}

if(!is.null(a$duong_scale)){
  df <- copy(XSRecency:::duong)
  df <- df[, days := days * as.numeric(a$duong_scale)]
  df <- df[, last.time := shift(days), by="id.key"]
  df <- df[, gap := days - last.time]
} else {
  df <- NULL
}

set.seed(a$seed)

if(!a$pt){
  sim <- simulate(n_sims=a$n_sims, n=a$n,
                  inc.function=inc.function,
                  infection.function=infection.function,
                  baseline_incidence=a$inc, prevalence=a$p, rho=a$rho,
                  phi.func=phi.func,
                  bigT=a$bigT, tau=a$tau, ext_FRR=a$ext_FRR,
                  ext_df=df,
                  max_FRR=a$max_FRR,
                  last_point=a$last_point)
} else {
  # THESE ARE THE PRIOR TEST SETTINGS

  ptest.dist <- function(n, u) runif(n, a$t_min, a$t_max)
  ptest.prob <- function(u) a$q

  if(a$mech2){
    GAMMA_PARMS <- c(1.57243557, 1.45286770, -0.02105187)
    ptest.dist2 <- function(n, u) u - rgengamma(n=n,
                                                mu=GAMMA_PARMS[1],
                                                sigma=GAMMA_PARMS[2],
                                                Q=GAMMA_PARMS[3])
  } else {
    ptest.dist2 <- NULL
  }

  if(!is.null(a$gamma)){
    t_noise <- function(t) min(0, t + rnorm(n=1, sd=a$gamma))
  } else {
    t_noise <- NULL
  }

  sim <- simulate.pt(n_sims=a$n_sims, n=a$n,
                     inc.function=inc.function,
                     infection.function=infection.function,
                     baseline_incidence=a$inc, prevalence=a$p, rho=a$rho,
                     phi.func=phi.func,
                     bigT=a$bigT, tau=a$tau, ext_FRR=a$ext_FRR,
                     ext_df=df,
                     max_FRR=a$max_FRR,
                     last_point=a$last_point,
                     # THESE ARE THE NEW ARGUMENTS --
                     ptest.dist=ptest.dist,
                     ptest.prob=ptest.prob,
                     t_range=c(a$t_min, a$t_max),
                     t_noise=t_noise,
                     d_misrep=a$eta,
                     q_misrep=a$nu,
                     p_misrep=a$xi,
                     ptest.dist2=ptest.dist2)
}

df <- do.call(cbind, sim) %>% data.table
df[, sim := .I]

as <- do.call(c, a)
for(i in 1:length(as)){
  df[, names(as[i]) := as[i]]
}
print(df)
print(df)

filename <- do.call(paste, a)
filename <- gsub(" ", "_", filename)
write.csv(df, file=paste0(out_dir, "results-", filename, ".csv"))
