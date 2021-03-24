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

fit.model <- function(df){

  # Fit GEE model for the gap times
  mod <- gee(gap ~ 1, id=id.key, data=df,
             family=poisson(link="log"), corstr="exchangeable")

  # Get only the first days of sampling
  first.day <- df[samp == 1]

  return(list(days=first.day$days,
              num.samples=first.day$num.samples,
              rate=exp(coef(mod))))
}

simulate.study <- function(days, num.samples, rate){

  # Sample the first day
  day1 <- sample(days, size=length(days), replace=TRUE)

  # Sample the number of samples
  nums <- sample(num.samples, size=length(days), replace=TRUE)

  # Get the rate for gap times
  # rate <- cbind(rep(1, length(days)), day1, nums) %*% coef(mod)

  # Get the gap-times for each
  gaps <- lapply(nums, function(x) rpois(n=x, lambda=rate))

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

df <- process.study(file=FILE)
num.samples <- df[, lapply(.SD, max), .SDcols="samp", by="id.key"]
samps <- num.samples$samp
pdf("duong-hist-data.pdf", width=8, height=4)
par(mfrow=c(1, 2))
hist(samps, main="Histogram of Number of Samples",
     xlab="Number of Samples Taken")
hist(df$days, main="Days Post Seroconversion Histogram",
     xlab="Days")
dev.off()

sims <- simulate.studies(file=FILE, nsims=100)

par(mfrow=c(2, 2))
hist(df$days)
hist(sims[[1]]$days)
hist(sims[[2]]$days)
hist(sims[[3]]$days)

# Fit GEE model for the gap times
df[, samp.5 := samp >= 5]
mod <- gee(gap ~ samp + samp.5 + samp.5*samp, id=id.key, data=df,
           family=poisson(link="log"), corstr="exchangeable")

dd <- data.table(samp=1:25)
dd[, samp.5 := samp >= 5]

int <- rep(1, 25)
samp <- 1:25
samp.5 <- samp >= 5

dd <- cbind(int, samp, samp.5, samp*samp.5)
preds <- exp(dd %*% coef(mod))
dd <- data.table(dd)
dd[, preds := preds]

ggplot(data=df) + geom_line(aes(x=samp, y=gap, group=id.key)) +
  geom_line(data=dd, aes(x=samp, y=preds), color='red')


# TODO: Put this into the function.
first.day <- df[samp==1]
days <- first.day$days
num.samples <- first.day$num.samples
# Sample the first day
day1 <- sample(days, size=length(days), replace=TRUE)

# Sample the number of samples
nums <- sample(num.samples, size=length(days), replace=TRUE)

# Get sample numbers
samp.num <- lapply(nums, function(x) 1:x)

get.samp.dmat <- function(samp) return(cbind(1, samp, samp>=5, (samp>=5)*samp))

dmats <- lapply(samp.num, get.samp.dmat)

rates <- lapply(dmats, function(x) exp(x %*% coef(mod))[, 1])

# Get the gap-times for each
gaps <- lapply(rates, function(x) rpois(n=length(x), lambda=x))

# Add on the start day
days.gap <- mapply(function(x, y) c(x, y), x=day1, y=gaps)

# Calculate sampling times
days.tot <- lapply(days.gap, cumsum)

# Get ids and put into data frame
ids <- 1:length(days)
ids <- rep(ids, nums)

hist(unlist(days.tot), breaks=50)




hist(df$days, main="Days Post Seroconversion Histogram",
     xlab="Days", breaks=50)
