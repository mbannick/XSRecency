rm(list=ls())

library(data.table)
library(magrittr)
library(stringr)
library(flexsurv)

# Load XSRecency
library(XSRecency)
library(simtools)

source("./sim-helpers.R")

if(is.parallel()){
  a <- commandArgs(trailingOnly=TRUE)
  OUTDIR <- a[1]
} else {
  OUTDIR <- "."
}

# TODO: Introduce TASK ID and param getter.

# Get the gamma parameters and baseline phi function
phi.func <- get.phi()

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
  if(!a$pt){
    stop("Piecewise constant-linear incidence function only to be used
         with prior testing simulations.")
  }
  infection.function <- function(...) infections.lincon(bigT=a$bigT, ...)
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

  ptest.dist <- function(u) runif(1, a$t_min, a$t_max)
  ptest.prob <- function(u) a$q

  if(a$mech2){
    GAMMA_PARMS <- c(1.57243557, 1.45286770, -0.02105187)
    ptest.dist2 <- function(u) u - rgengamma(n=1,
                                             mu=GAMMA_PARMS[1],
                                             sigma=GAMMA_PARMS[2],
                                             Q=GAMMA_PARMS[3])
  } else {
    ptest.dist2 <- NULL
  }

  if(!is.null(a$gamma)){
    t_noise <- function(t) max(0, t + rnorm(n=1, sd=a$gamma))
  } else {
    t_noise <- NULL
  }

  if(!is.null(a$t_min_exclude)){
    t_min <- a$t_min_exclude
  } else {
    t_min <- 0
  }
  start <- Sys.time()
  sim <- simulate.pt(n_sims=a$n_sims, n=a$n,
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
                     t_max=100,
                     t_min=t_min,
                     t_noise=t_noise,
                     d_misrep=a$eta,
                     q_misrep=a$nu,
                     p_misrep=a$xi,
                     ptest.dist2=ptest.dist2,
                     exclude_pt_bigT=a$exclude_pt_bigT,
                     t_total_exclude=a$t_total_exclude)
  end <- Sys.time()
  print(end - start)
}

df <- do.call(cbind, sim) %>% data.table
df[, sim := .I]

as <- do.call(c, a)
for(i in 1:length(as)){
  if((names(as[i])) %in% colnames(df)) next
  df[, names(as[i]) := as[i]]
}
print(df)

filename <- do.call(paste, a)
filename <- gsub(" ", "_", filename)
write.csv(df, file=paste0(out_dir, "results-", filename, ".csv"))
