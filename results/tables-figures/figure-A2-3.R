rm(list=ls())

library(data.table)
library(ggplot2)
library(geepack)
library(lme4)
library(gridExtra)
library(magrittr)
library(XSRecency)

# This data was obtained from:
# https://doi.org/10.1371/journal.pone.0114947.s001
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
sims <- simulate.studies(nsims=3)

pdf("duong-hist-fig.pdf", width=8, height=5)
par(mfrow=c(1, 2))
hist(XSRecency:::duong$days/365.25, main="Years Post Seroconversion \n(Observed)",
     xlab="Years", breaks=50, ylim=c(0, 205))
hist(sims[[1]]$durations, main="Years Post Seroconversion \n(One Simulation)",
     xlab="Years", breaks=50, col='lightblue',
     ylim=c(0, 205))
dev.off()

# EXPLANATION OF METHOD ---------------

df <- copy(duong)
# df[, days := days * 10/8]
df[, last.time := shift(days), by="id.key"]
df[, gap := days - last.time]

ggplot(data=df, aes(y=gap, x=num.samples)) + geom_jitter() + theme_bw()

df[, samp.5 := samp >= 5]
mod <- geese(gap ~ samp + samp.5 + samp.5*samp, id=id.key, data=df,
             family=poisson(link="log"), corstr="exchangeable")

samp.num <- 1:25
get.samp.dmat <- function(samp) return(cbind(1, samp, samp>=5, (samp>=5)*samp))
dmats <- lapply(samp.num, get.samp.dmat)

# Compute the predicted rates
rates <- lapply(dmats, function(x) exp(x %*% mod$beta)[, 1]) %>% unlist
df2 <- data.frame(samp=samp.num, gaps=rates)

df1 <- df[, list(id.key, first.samp)] %>% unique()
df3 <- df[, list(id.key, num.samples)] %>% unique()

pdf("sample-times.pdf", height=5, width=10)
p1 <- ggplot(data=df, aes(y=gap, x=samp)) + geom_jitter(alpha=0.5) +
  geom_line(data=df2, aes(y=gaps, x=samp), color='red') +
  theme_minimal() + xlab("Sample Number") + ylab("Gap Time Between Samples (Days)")
p2 <- ggplot(data=df1) +
  geom_histogram(aes(first.samp), color="black", fill="grey", bins=10) +
  xlab("Infection Duration (Days) at First Sample") + ylab("Frequency") +
  theme_minimal()
p3 <- ggplot(data=df3) +
  geom_histogram(aes(num.samples), color="black", fill="grey", bins=10) +
  xlab("Total Number of Samples") + ylab("Frequency") +
  theme_minimal()
grid.arrange(p2, p3, p1, nrow=1)
dev.off()
