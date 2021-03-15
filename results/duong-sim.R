rm(list=ls())

library(data.table)
library(ggplot2)
library(gee)
library(lme4)
library(gridExtra)
library(magrittr)

# --------------------------------------------------------------------
# DATA CLEANING ------------------------------------------------------
# --------------------------------------------------------------------

df <- fread("~/OneDrive/Documents/2020_2021/RA FEI GAO/duong2015.csv")
df <- df[!is.na(days)]

ids <- df[, .(id1)]
ids <- unique(ids)
ids[, id.key := .I]

df <- merge(df, ids, by=c("id1"))
setorder(df, id.key, days)
df <- df[, .(id.key, days)]
df[, samp := 1:.N, by="id.key"]
df[days == 0, days := 0.5]

hist(df$days, main="Days Post Sero-Conversion Histogram",
     xlab="Days")

# --------------------------------------------------------------------
# TRAJECTORIES OF SAMPLING TIMES WITHIN A PANEL ----------------------
# --------------------------------------------------------------------

# Plots of the days since seroconversion of sampling over the sampling id
p1 <- ggplot(data=df) + geom_line(aes(x=samp, y=days, group=id.key)) +
  ggtitle("Trajectories of Days Post\nSero-Conversion by Sample Number") +
  labs(x="Sample Number", y="Days Post Sero-Conversion")
p2 <- ggplot(data=df) + geom_line(aes(x=samp, y=log(days), group=id.key)) +
  ggtitle("Trajectories of Log Days Post\nSero-Conversion by Sample Number") +
  labs(x="Sample Number", y="Log Days Post Sero-Conversion")
grid.arrange(p1, p2, nrow=1)

# Make individual trajectories of how many
# days post seroconversion it took they were sampled at
df[, samp2 := samp**2]
df[, samp3 := samp**3]
mod <- lmer(log(days) ~ samp + samp2 + samp3 + (1|id.key), data=df)
# mod <- lmer(log(days) ~ samp + (1|id.key), data=df)
df[, pred := predict(mod)]

pdf("duong-trajectories-data.pdf", width=6, height=6)
ggplot(data=df) + geom_line(aes(x=samp, y=days, group=id.key)) +
  # geom_line(aes(x=samp, y=exp(pred), group=id.key), color='red') +
  ggtitle("Individual Trajectories") +
  labs(x="Sample Number", y="Days Post Seroconversion")
dev.off()

# --------------------------------------------------------------------
# NUMBER OF SAMPLES TAKEN PER PANEL ----------------------------------
# --------------------------------------------------------------------

num.samples <- df[, lapply(.SD, max), .SDcols="samp", by="id.key"]$samp
pdf("duong-hist-data.pdf", width=8, height=4)
par(mfrow=c(1, 2))
hist(num.samples, main="Histogram of Number of Samples",
     xlab="Number of Samples Taken")
hist(df$days, main="Days Post Seroconversion Histogram",
     xlab="Days")
dev.off()

# --------------------------------------------------------------------
# DATA GENERATING MECHANISM ------------------------------------------
# --------------------------------------------------------------------

ids <- unique(df$id.key)

make.trajectory <- function(id, n){
  cos <- coef(mod)$id.key[id, ] %>% t
  samp.id <- 1:n
  variables <- cbind(rep(1, n), samp.id, samp.id**2, samp.id**3)
  # variables <- cbind(rep(1, n), samp.id)
  days <- exp(variables %*% cos) %>% round
  return(days)
}

simulate.one <- function(){

  total.samples <- sample(x=num.samples, size=length(ids), replace=TRUE)
  traj <- mapply(make.trajectory, id=ids, n=total.samples, SIMPLIFY=T)

  id.sims <- rep(ids, total.samples)
  times <- unlist(traj)

  new.df <- data.table(id=id.sims, days=times)
}

set.seed(10)

ps <- list()

pdf("duong-hist.pdf", width=10, height=10)
par(mfrow=c(4, 2))
hist(num.samples, main="DATA: Histogram of \nNumber of Samples",
     xlab="Number of Samples Taken", xlim=c(0, 30), breaks=30,
     col='pink')
hist(df$days, main="DATA: Days Post \nSeroconversion Histogram",
     xlab="Days Post Seroconversion", xlim=c(0, 1e4), breaks=seq(0, 1e4, by=50),
     col='pink')
p0 <- ggplot(data=df) + geom_line(aes(x=samp, y=days, group=id.key), color="#ff59c8") +
  ggtitle("DATA: Individual Trajectories") + ylim(c(0, 10000)) +
  ylab("Days Post Seroconversion") + xlab("Sample Number")
ps[[1]] <- p0

for(i in 1:3){

  sim <- simulate.one()
  sim[, samp := 1:.N, by="id"]

  sim.samples <- sim[, .N, by="id"]
  hist(sim.samples$N, main="SIM: Histogram of \nNumber of Samples",
       xlab="Number of Samples Taken", xlim=c(0, 30), breaks=30)
  hist(sim$days, main="SIM: Days Post \nSeroconversion Histogram",
       xlab="Days Post Seroconversion", xlim=c(0, 1e4), breaks=seq(0, 1e4, by=50))

  p.i <- ggplot(data=sim) + geom_line(aes(x=samp, y=days, group=id)) +
    ggtitle("SIM: Individual Trajectories") + ylim(c(0, 10000)) +
    ylab("Days Post Seroconversion") + xlab("Sample Number")
  ps[[i+1]] <- p.i
}
dev.off()

pdf("duong-trajectories.pdf", width=6, height=6)
grid.arrange(grobs=ps)
dev.off()
