rm(list=ls())

library(data.table)
library(ggplot2)
library(gee)
library(lme4)
library(gridExtra)
library(magrittr)
library(XSRecency)

FILE <- "data-raw/duong2015.csv"

# --------------------------------------------------------------------
# DATA GENERATING MECHANISM ------------------------------------------
# --------------------------------------------------------------------

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
  df <- df[, .(id.key, days)]
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

# Raw Data Analysis --------------------------------------------
# Save for internal package functions
# duong <- process.study(file=FILE)
# usethis::use_data(duong, internal=TRUE, overwrite=TRUE)

# Simulation ---------------------------------------------------
set.seed(365)
sims <- simulate.studies(file=FILE, nsims=3)

pdf("duong-hist-fig.pdf", width=8, height=5)
par(mfrow=c(1, 2))
hist(XSRecency:::duong$days, main="Days Post Seroconversion \n(Observed)",
     xlab="Days", breaks=50, ylim=c(0, 205))
hist(sims[[1]]$days, main="Days Post Seroconversion \n(One Simulation)",
     xlab="Days", breaks=50, col='lightblue',
     ylim=c(0, 205))
dev.off()
