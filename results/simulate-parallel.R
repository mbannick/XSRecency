rm(list=ls())

args <- commandArgs()
print(args)

library(data.table)
library(magrittr)
library(R.utils)
library(stringr)

setwd("~/repos/XSRecency/")
source("./results/sim-helpers.R")
source("./R/phi-functions.R")

# Get command-line arguments
a <- commandArgs(trailingOnly=TRUE, asValues=TRUE,
                    defaults=list(
                      seed=1,
                      n_sims=10,
                      n=5000,
                      p=0.29,
                      inc=0.032,
                      window=248,
                      shadow=306,
                      itype="constant",
                      rho=NULL,
                      tau=12,
                      bigT=2,
                      phi_frr=0.02,
                      phi_tfrr=NULL,
                      phi_norm_mu=7,
                      phi_norm_sd=1,
                      phi_norm_div=8,
                      out_dir=".",
                      ext_FRR=FALSE,
                      duong_scale=NULL,
                      max_FRR=NULL,
                      last_point=FALSE
                    ))

# Capture date in the out directory
date <- format(Sys.time(), "%d-%m-%y-%H")
out_dir <- paste0(a$out_dir, "/", date, "/")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# a[[1]] <- NULL
a$out_dir <- NULL

print(a)

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
phit.pnorm <- function(t) phit.const(t) + pnorm(t, mean=6.5, sd=1) / 8

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

sim <- simulate(n_sims=a$n_sims, n=a$n,
                inc.function=inc.function,
                infection.function=infection.function,
                baseline_incidence=a$inc, prevalence=a$p, rho=a$rho,
                phi.func=phi.func,
                bigT=a$bigT, tau=a$tau, ext_FRR=a$ext_FRR,
                ext_df=df,
                max_FRR=a$max_FRR,
                last_point=a$last_point)

df <- do.call(cbind, sim) %>% data.table
df[, sim := .I]

as <- do.call(c, a)
for(i in 1:length(as)){
  df[, names(as[i]) := as[i]]
}

filename <- do.call(paste, a)
filename <- gsub(" ", "_", filename)
write.csv(df, file=paste0(out_dir, "results-", filename, ".csv"))
