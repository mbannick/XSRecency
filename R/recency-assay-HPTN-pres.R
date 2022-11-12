setwd("~/OneDrive/Documents/2020-2021/RA FEI GAO/HPTN-faculty-Oct-3/")

# PHI FIGURE ----------------------------------

pdf('assay.pdf',height=6,width=6)
par(mfrow=c(1,1))

mdri <- 101
shadow <- 194

t = seq(0, 12, 1e-3)
params <- get.gamma.params(window=mdri/365.25, shadow=shadow/365.25)

phit <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

i <- 1

if(i == 1) ttime <- 2
if(i == 1) tval <- phit(2)
if(i == 1) name <- "Assay 1"

if(i == 2) tval <- 0.02
if(i == 2) ttime <- uniroot(function(t) phit(t) - tval, interval=c(0, 12))$root
if(i == 2) name <- "Assay 2"

phit.const <- function(t) phit(t)*(t <= ttime) + tval*(t > ttime)
# TRY THIS FUNCTION
# plot(t, phit(t), col='blue', type='l', ylab=expression(phi(t)), xlab=expression(t))

###### MDRI at 2 year 142/365, FRR 1.5%: 1.5% + Gamma distribution with mean 142/365.25-1.5%
plot(t, phit.const(t), col='red', type='l', ylab=expression(phi(t)), xlab=expression(t))
abline(h=tval, lty='dashed')

if(i == 2){
  pgon <- c(seq(2, ttime, 0.01))
  pgon.phi <- 1-pgamma(pgon, shape = params[1], rate = params[2])
  polygon(x=c(pgon, rev(pgon)), y=c(rep(tval, length(pgon)), rev(pgon.phi)), col='lightgrey')
}
abline(v=2, lty='dashed')

dev.off()

pdf('assay.pdf',height=6,width=10)
par(mfrow=c(1,2))

for(i in 1:2){

  if(i == 1){
    mdri <- 101
    shadow <- 194
  } else {
    mdri <- 248
    shadow <- 306
  }

  t = seq(0, 12, 1e-3)
  params <- get.gamma.params(window=mdri/365.25, shadow=shadow/365.25)

  phit <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])

  if(i == 1) ttime <- 2
  if(i == 1) tval <- phit(2)
  if(i == 1) name <- "Assay 1"

  if(i == 2) tval <- 0.02
  if(i == 2) ttime <- uniroot(function(t) phit(t) - tval, interval=c(0, 12))$root
  if(i == 2) name <- "Assay 2"

  phit.const <- function(t) phit(t)*(t <= ttime) + tval*(t > ttime)
  # TRY THIS FUNCTION
  plot(t, phit(t), col='blue', type='l', ylab=expression(phi(t)), xlab=expression(t))
  abline(h=tval, lty='dashed')

  ###### MDRI at 2 year 142/365, FRR 1.5%: 1.5% + Gamma distribution with mean 142/365.25-1.5%
  lines(t, phit.const(t), col='red')

  if(i == 2){
    pgon <- c(seq(2, ttime, 0.01))
    pgon.phi <- 1-pgamma(pgon, shape = params[1], rate = params[2])
    polygon(x=c(pgon, rev(pgon)), y=c(rep(tval, length(pgon)), rev(pgon.phi)), col='lightgrey')
  }

  legend('topright', paste0(name, c("A", "B")),
         lty=rep(1, 4),col=c("blue", "red"), cex=0.75)
  abline(v=2, lty='dashed')
}

dev.off()

library(data.table)
library(ggplot2)

df <- fread("data-raw/duong2015-with-lag.csv")
ids <- df[, .(id1)]
ids <- unique(ids)
ids[, id.key := .I]

df <- merge(df, ids, by=c("id1"))
setorder(df, id.key, days)
df <- df[, .(id.key, days, LAg)]
df[, samp := 1:.N, by="id.key"]

pdf("duong.pdf", height=5, width=8)
ggplot(data=df) +
  geom_line(aes(x=days, y=LAg, group=id.key), alpha=0.5) +
  geom_point(aes(x=days, y=LAg), alpha=1) +
  labs(x="Days Post Seroconversion") +
  # geom_hline(yintercept=1.5) +
  geom_vline(xintercept=365.25, linetype='dashed')
dev.off()

df$recent <- df$LAg <= 1.5

pdf("duong-phi.pdf", height=5, width=8)
par(mfrow=c(1, 2))

t <- seq(0, 3000, 1)
y <- as.numeric(t <= 365.25)
plot(y ~ t, type='l', ylab="Test Recent Probability",
     xlab="Days Post Seroconversion")
abline(v=365.25, lty='dashed')

mod <- with(df, fit.cubic(recent=recent, durations=days, id=id.key))
d <- seq(0, max(df$days, na.rm=T), by=1)
phihat <- phi.hat(d, mod)
datt <- data.frame(d=d, phi=c(phihat))

library(scales)
plot(as.numeric(df$recent) ~ df$days, pch=16,
     col=alpha("grey", 0.7),
     ylab="Test Recent Probability",
     xlab="Days Post Seroconversion")
with(datt, lines(phihat ~ d))
abline(v=365.25, lty='dashed')
dev.off()

ggplot(data=datt, mapping=aes(x=d, y=phihat)) +
  geom_point(data=df[!is.na(LAg)], aes(x=days, y=as.numeric(recent)),
             alpha=0.2) +
  geom_line(color="red") +
  labs(x="Days Post Seroconversion",
       y="Test Recent Probability") +
  geom_vline(xintercept=365.25, linetype='dashed')

integrate.phi(mod, minT=0, maxT=1)
mean(df[days > 365.25, recent], na.rm=T)

df[recent == TRUE & days <= 365.25, col := "True Recent"]
df[recent == TRUE & days > 365.25, col := "False Recent"]
df[recent == FALSE & days <= 365.25, col := "False Long"]
df[recent == FALSE & days > 365.25, col := "True Long"]

pdf("duong-quad-1.pdf", height=5, width=8)
ggplot(data=df[col %in% c("True Recent", "True Long") &
                 !is.na(LAg)]) +
  geom_point(aes(x=days, y=LAg, color=col), alpha=0.7) +
  labs(x="Days Post Seroconversion") +
  geom_vline(xintercept=365.25, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  scale_color_manual(values=c("#ff3b3b", "#008be8")) +
  theme(legend.title=element_blank(),
        legend.position="top") +
  ylim(c(min(df$LAg, na.rm=T), max(df$LAg, na.rm=T)))
dev.off()

pdf("duong-quad-2.pdf", height=5, width=8)
ggplot(data=df[!is.na(LAg)]) +
  geom_point(aes(x=days, y=LAg, color=col), alpha=0.7) +
  labs(x="Days Post Seroconversion") +
  geom_vline(xintercept=365.25, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  scale_color_manual(values=c("#00e808", "#e89700",
                              "#ff3b3b", "#008be8")) +
  theme(legend.title=element_blank(),
        legend.position="top") +
  ylim(c(min(df$LAg, na.rm=T), max(df$LAg, na.rm=T)))
dev.off()

df$pt <- rbinom(nrow(df), prob=0.1, size=1) %>% as.factor()
df$ly <- rbinom(nrow(df), prob=0.8, size=1) %>% as.factor()
df$use <- ((df$pt == 1) & (df$ly == 1) & (df$days <= 365.25)) %>% as.factor()

df[, pt_char := factor(pt, labels=c("No", "Yes"))]
df[, use_char := factor(use, labels=c("No", "Yes"))]
df[, recent_new := (recent | as.logical(use))]

pdf("duong-quad-3.pdf", height=5, width=8)
ggplot() +
  geom_point(data=df[!is.na(LAg)],
             aes(x=days, y=LAg,
                 alpha=pt_char,
                 shape=pt_char,
                 size=pt_char), color="black") +
  labs(x="Days Post Seroconversion",
       alpha="Has Negative Prior Test",
       shape="Has Negative Prior Test",
       size="Has Negative Prior Test") +
  geom_vline(xintercept=365.25, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  scale_alpha_manual(values=c(0.1, 1.0)) +
  scale_shape_manual(values=c("circle", "diamond plus")) +
  scale_size_manual(values=c(1, 3)) +
  theme(legend.position="top") +
  ylim(c(min(df$LAg, na.rm=T), max(df$LAg, na.rm=T)))
dev.off()

pdf("duong-quad-4.pdf", height=5, width=8)
ggplot() +
  geom_point(data=df[!is.na(LAg)],
             aes(x=days, y=LAg,
                 alpha=use_char,
                 shape=use_char,
                 size=use_char), color="black") +
  labs(x="Days Post Seroconversion",
       alpha="Has RECENT Negative Prior Test",
       shape="Has RECENT Negative Prior Test",
       size="Has RECENT Negative Prior Test") +
  geom_vline(xintercept=365.25, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  scale_alpha_manual(values=c(0.1, 1.0)) +
  scale_shape_manual(values=c("circle", "diamond plus")) +
  scale_size_manual(values=c(1, 3)) +
  theme(legend.position="top") +
  ylim(c(min(df$LAg, na.rm=T), max(df$LAg, na.rm=T)))
dev.off()

df[recent_new == TRUE & days <= 365.25, col := "True Recent"]

pdf("duong-quad-5.pdf", height=5, width=8)
ggplot(data=df[!is.na(LAg)]) +
  geom_point(aes(x=days, y=LAg, color=col), alpha=0.7) +
  labs(x="Days Post Seroconversion") +
  geom_vline(xintercept=365.25, linetype='dashed') +
  geom_hline(yintercept=1.5, linetype='dashed') +
  scale_color_manual(values=c("#00e808", "#e89700",
                              "#ff3b3b", "#008be8")) +
  theme(legend.title=element_blank(),
        legend.position="top") +
  ylim(c(min(df$LAg, na.rm=T), max(df$LAg, na.rm=T)))
dev.off()

mdri <- 248
shadow <- 306
params <- get.gamma.params(window=mdri/365.25, shadow=shadow/365.25)
phit <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
ts <- seq(0, 8, 1e-3)

pdf("~/OneDrive/Documents/2020-2021/RA FEI GAO/HPTN-faculty-Oct-3/phi.pdf",
    height=5, width=5)
par(mfrow=c(1, 1))
plot(phit(ts) ~ ts, type='l', ylab=expression(phi(t)), xlab=expression(t))
dev.off()

pdf("~/OneDrive/Documents/2020-2021/RA FEI GAO/HPTN-faculty-Oct-3/meanwindow.pdf",
    height=5, width=5)
par(mfrow=c(1, 1))
plot(phit(ts) ~ ts, type='l', ylab=expression(phi(t)), xlab=expression(t))
pgon <- c(seq(0, 8, 0.01))
pgon.phi <- 1-pgamma(pgon, shape = params[1], rate = params[2])
polygon(x=c(pgon, rev(pgon)), y=c(rep(0, length(pgon)), rev(pgon.phi)), col='#5b92f0')
dev.off()

pdf("~/OneDrive/Documents/2020-2021/RA FEI GAO/HPTN-faculty-Oct-3/mdri.pdf",
    height=5, width=5)
par(mfrow=c(1, 1))
plot(phit(ts) ~ ts, type='l', ylab=expression(phi(t)), xlab=expression(t))
pgon <- c(seq(0, 2, 0.01))
pgon2 <- c(seq(2, 8, 0.01))
pgon.phi <- 1-pgamma(pgon, shape = params[1], rate = params[2])
pgon.phi2 <- 1-pgamma(pgon2, shape=params[1], rate=params[2])
polygon(x=c(pgon, rev(pgon)), y=c(rep(0, length(pgon)), rev(pgon.phi)), col='#5b92f0')
polygon(x=c(pgon2, rev(pgon2)), y=c(rep(0, length(pgon2)), rev(pgon.phi2)), col='#f05b5b')
dev.off()
