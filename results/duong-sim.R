rm(list=ls())

library(data.table)
library(ggplot2)
library(gee)
library(lme4)
library(gridExtra)
library(magrittr)

FILE <- "~/OneDrive/Documents/2020_2021/RA FEI GAO/duong2015.csv"

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

fit.model <- function(df, knot=5){

  # Indicator variable for sample 5+
  df[, samp.5 := samp >= knot]

  # Fit GEE model for the gap times
  mod <- gee(gap ~ samp + samp.5 + samp.5*samp, id=id.key, data=df,
             family=poisson(link="log"), corstr="exchangeable")

  # Get only the first days of sampling
  first.day <- df[samp == 1]

  return(list(days=first.day$days,
              num.samples=first.day$num.samples,
              rate=coef(mod)))
}

simulate.study <- function(days, num.samples, coefs, knot=5){
  # Sample the first day
  day1 <- sample(days, size=length(days), replace=TRUE)

  # Sample the number of samples
  nums <- sample(num.samples, size=length(days), replace=TRUE)

  # Get a vector of the sample numbers, starting from the second one
  samp.num <- lapply(nums, function(x) 2:x)

  # Get a sample design matrix
  get.samp.dmat <- function(samp) return(cbind(1, samp, samp>=knot, (samp>=knot)*samp))

  dmats <- lapply(samp.num, get.samp.dmat)

  rates <- lapply(dmats, function(x) exp(x %*% coefs)[, 1])

  # Get the gap-times for each
  gaps <- lapply(rates, function(x) rpois(n=length(x), lambda=x))

  # Add on the start day
  days.gap <- mapply(function(x, y) c(x, y), x=day1, y=gaps)

  # Calculate sampling times
  days.tot <- lapply(days.gap, cumsum)

  # Get ids and put into data frame
  ids <- 1:length(days)
  ids <- rep(ids, nums)

  df <- data.table(id.key=ids,
                   days=unlist(days.tot))
  return(df)
}

simulate.studies <- function(file, nsims){
  df <- process.study(file=file)
  mod <- fit.model(df=df)
  sims <- replicate(nsims, simulate.study(mod$days, mod$num.samples, mod$rate),
                    simplify=FALSE)
  return(sims)
}

# --------------------------------------------------------------------
# NUMBER OF SAMPLES TAKEN PER PANEL ----------------------------------
# --------------------------------------------------------------------

# Raw Data
df <- process.study(file=FILE)
num.samples <- df[, lapply(.SD, max), .SDcols="samp", by="id.key"]
samps <- num.samples$samp

set.seed(365)
# Simulation
sims <- simulate.studies(file=FILE, nsims=3)

pdf("duong-hist-fig.pdf", width=8, height=5)
par(mfrow=c(1, 2))
hist(df$days, main="Days Post Seroconversion \n(Observed)",
     xlab="Days", breaks=50, ylim=c(0, 205))
hist(sims[[1]]$days, main="Days Post Seroconversion \n(One Simulation)",
     xlab="Days", breaks=50, col='lightblue',
     ylim=c(0, 205))
dev.off()
