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
                      window=71,
                      shadow=80,
                      itype="constant",
                      rho=NULL,
                      tau=12,
                      bigT=2,
                      phi_frr=NULL,
                      phi_tfrr=NULL,
                      phi_norm_mu=NULL,
                      phi_norm_sd=NULL,
                      phi_norm_div=NULL,
                      out_dir="."
                    ))

print(a)

# Capture date in the out directory
date <- format(Sys.time(), "%d-%m-%y-%h")
out_dir <- paste0(a$out_dir, "/", date, "/")
dir.create(out_dir, showWarnings=FALSE)

# a[[1]] <- NULL
a$out_dir <- NULL

print(a)

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
if(is.null(a$rho)) rho <- NA

# Get the gamma parameters and baseline phi function
params <- get.gamma.params(window=a$window/356.25, shadow=a$shadow/365.25)

# Set up each type of phi function, will be overwritten
phi.none <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
phi.const <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
phi.dnorm <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

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

if(a$itype == "constant"){
  inc.function <- c.incidence
  infection.function <- c.infections
} else if(a$itype == "linear"){
  inc.function <- l.incidence
  infection.function <- l.infections
} else if(a$itype == "exponential"){
  inc.function <- e.incidence
  infection.function <- e.infections
} else {
  stop("Unknown incidence function.")
}

set.seed(a$seed)

sim <- simulate(n_sims=a$n_sims, n=a$n,
                inc.function=inc.function,
                infection.function=infection.function,
                baseline_incidence=a$inc, prevalence=a$p, rho=a$rho,
                phi.func=phi.func,
                bigT=a$bigT, tau=a$tau)

df <- do.call(cbind, sim) %>% data.table
df[, sim := .I]

as <- do.call(c, a)
for(i in 1:length(as)){
  df[, names(as[i]) := as[i]]
}

filename <- do.call(paste, a)
filename <- gsub(" ", "_", filename)
write.csv(df, file=paste0(out_dir, "results-", filename, ".csv"))
