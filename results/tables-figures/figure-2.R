rm(list=ls())

library(data.table)
library(ggplot2)
library(geepack)
library(lme4)
library(gridExtra)
library(magrittr)

setwd("~/repos/XSRecency/")
library(XSRecency)

# This data was obtained from:
# https://doi.org/10.1371/journal.pone.0114947.s001
FILE <- "data-raw/duong2015-with-lag.csv"

# PROCESS DUONG 2015 DATA -------------------------
process.study <- function(file){
  # Read in data
  df <- fread(file)

  # Data pre-processing
  df <- df[!is.na(days)]

  ids <- df[, .(id1)]
  ids <- unique(ids)
  ids[, id.key := .I]

  df <- merge(df, ids, by=c("id1"))
  setorder(df, id.key, days)
  df <- df[, .(id.key, days, LAg, type)]
  df[, samp := 1:.N, by="id.key"]

  # Calculate gap days between samples
  df[, last.time := shift(days), by="id.key"]
  df[, gap := days - last.time]

  # Calculate total number of samples
  df[, num.samples := .N, by="id.key"]

  # Calculate first day
  df[, first.samp := lapply(.SD, min), .SDcols="days", by="id.key"]

  return(df)
}

duong <- process.study(FILE)
frr <- fread("data-raw/duong2015-frr.csv")
BIGT <- 1
TAU <- 8

study.data <- duong[, .(id.key, days, LAg)]
study.data[, recent := as.integer(LAg <= 1.5)]
study.data[, durations := days / 365.25]
setnames(study.data, c("id.key"), c("id"))

# Calculate window period and MDRI based on T^* and \tau
windows <- assay.properties.est(study=study.data, bigT=BIGT, tau=TAU)

study.data.frr <- study.data[days > BIGT*365.25 & !is.na(recent)]

# Calculate false recency rate
frr[, recent := LAg <= 1.5]
frr <- frr[!is.na(LAg)] # ID number 6085 was missing LAg
beta <- sum(frr$recent) / nrow(frr) # this is really high
beta_var <- beta * (1 - beta) / nrow(frr)

# Data from the Negedu-Momoh et al. 2021 paper
N <- 8306
N_pos <- 394
N_neg <- N - N_pos
N_rec <- 19

snapshot <- get.snapshot(n=N, n_r=N_rec, n_n=N_neg, n_p=N_pos,
                         mu=windows$mu_est, mu_var=windows$mu_var)
adjusted <- get.adjusted(n=N, n_r=N_rec, n_n=N_neg, n_p=N_pos,
                         omega=windows$omega_est, omega_var=windows$omega_var,
                         beta=beta, beta_var=beta_var, big_T=BIGT)

snap <- snapshot$est + c(-1, 0, 1) * qnorm(0.975) * snapshot$var**0.5
adj <- adjusted$est + c(-1, 0, 1) * qnorm(0.975) * adjusted$var**0.5

print(snap)
print(adj)
