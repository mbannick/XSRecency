rm(list=ls())
library(pbapply)
library(pbmcapply)
library(data.table)
source("./results/sim-helpers.R")
source("./R/phi-functions.R")

set.seed(100)

NSIM <- 5000
PREVALENCE <- 0.29
INCIDENCE <- 0.032
N <- 5000

INCFUNC <- incidence.con
INFFUNC <- infections.con

WINDOW <- 101
SHADOW <- 194

Q <- 1.0
TMIN <- 0
TMAX <- 4

params <- get.gamma.params(window=WINDOW/365.25, shadow=SHADOW/365.25)
phi <- function(t) 1-pgamma(t, shape = params[1], rate = params[2])
tval <- phi(2)
phi.const <- function(t, ...) phi(t)*(t <= 2) + tval*(t > 2)

BETA <- tval

PHI <- phi.const
OMEGA <- integrate(PHI, lower=0, upper=2)$value

# one simulation

one.study <- function(i){
  Np <- rbinom(n=1, size=N, prob=PREVALENCE)
  Nn <- N - Np

  Ei <- runif(Np)
  Ui <- INFFUNC(e=Ei, t=0, p=PREVALENCE, lambda_0=INCIDENCE, rho=0)

  Ri <- rbinom(n=Np, size=1, prob=PHI(-Ui))

  Qi <- rbinom(Np, size=1, prob=Q)
  Ti <- runif(Np, min=TMIN, max=TMAX)
  Ai <- Qi & (Ti <= 2)
  Bi <- Qi & (Ti > 2)
  Oi <- sapply(Ti[Ai], function(x) integrate(PHI, lower=0, upper=x)$value, simplify=T)

  Di <- as.integer(-Ui >= Ti)
  Mi <- Ri * (1 - Bi * Di) + (1 - Ri) * Ai * (1 - Di)

  pr <- sum(Ri) / Np
  pr_star <- sum(Mi) / Np

  beta.est <- BETA
  # NBETA <- round(runif(n=1, min=1000, max=2000))
  # beta.est <- mean(rbinom(n=NBETA, size=1, prob=BETA))

  # Term 1
  w41 <- sum(Ai * Ti)
  w42 <- sum(Oi)
  w4 <- w41 - w42
  w5 <- beta.est * sum(Bi * Ti)

  w1 <- sum(Mi) - beta.est * sum(1-Bi)

  pA <- mean(Ai)
  pB <- mean(Bi)
  p <- Np / N
  muTA <- mean(Ti[Ai])
  varTA <- var(Ti[Ai])
  muTB <- mean(Ti[Bi])
  varTB <- var(Ti[Bi])
  omegaTA <- mean(Oi)
  omegaTAvar <- var(Oi)
  omegaTAstar <- mean(Ti[Ai] * Oi)

  cov45 <- N * p**2 * beta.est * pB * pA * muTB * (omegaTA - muTA)

  var41 <- N * p * pA * (muTA**2 *    (1-p*pA) + varTA)
  var42 <- N * p * pA * (omegaTA**2 * (1-p*pA) + omegaTAvar)
  var43 <- N * p * pA * (omegaTAstar - p * pA * muTA * omegaTA)

  var4 <- N * p * pA * (
    (varTA + omegaTAvar) + (muTA**2 + omegaTA**2) -
      p*pA * (muTA - omegaTA)**2 - 2 * omegaTAstar
  )
  var4 <- N * p * p_A * (
    (var_TA + omega_TA_var) + (mu_TA**2 + omega_TA**2) -
      p * p_A * (mu_TA - omega_TA)**2 +
      r_TA + r_TAprime * (n * p - p) - 2 * omega_TAstar
  )

  cov14 <- N * p * pA * (
    (muTA - omegaTA) * (pr - p * pr_star - beta.est * (1 - p + p*pB)) +
    INCIDENCE * (1-p) / p * (muTA**2 + varTA + omegaTA**2 + omegaTAvar - 2 * omegaTAstar)
  )

  w11 <- sum(Mi)
  w12 <- beta.est * sum(1-Bi)

  cov14.1 <- N * p * (mean(Mi * Ai * Ti)  - p * pr_star * pA * muTA)
  cov14.2 <- N * p * (sum(Mi[Ai] * Oi)/Np - p * pr_star * pA * omegaTA)
  cov14.3 <- N * p * beta.est * pA * muTA * (1-p * (1-pB))
  cov14.4 <- N * p * beta.est * pA * omegaTA * (1 - p * (1-pB))

  cov14.intermediate <- cov14.1 - cov14.2 - cov14.3 + cov14.4

  # Comparing mean(Mi * Ai * Ti)
  EMiAiTi <- mean(Mi * Ai * Ti)
  EMAT.1 <- pr
  # EMAT.1 <- beta.est + (OMEGA - beta.est * 2) * INCIDENCE * (1-p) / p
  EMAT.1 <- EMAT.1 * pA
  EMAT.1 <- EMAT.1 * (muTA)
  EMAT.2 <- INCIDENCE * (1-p) / p * pA * (muTA**2 + varTA - omegaTAstar)
  # We should have approximately
  # EMiAiTi ~~ EMAT.1 + EMAT.2
  #
  EMAO.1 <- pr
  # EMAO.1 <- beta.est + (OMEGA - beta.est * 2) * INCIDENCE * (1-p) / p
  EMAO.1 <- EMAO.1 * pA
  EMAO.1 <- EMAO.1 * (omegaTA)
  EMAO.2 <- INCIDENCE * (1-p) / p * pA * (omegaTAstar - omegaTA**2 - omegaTAvar)

  cov14.1.2 <- N * p * ((EMAT.1 + EMAT.2) - p * pr_star * pA * muTA)
  cov14.2.2 <- N * p * ((EMAO.1 + EMAO.2) - p * pr_star * pA * omegaTA)

  cov14.intermediate.2 <- cov14.1.2 - cov14.2.2 - cov14.3 + cov14.4

  EMiAiOi <- sum(Mi[Ai] * Oi)/Np

  EMi.1 <- mean(Ri)
  EMiAi.1 <- mean(Ri * Ai)
  EMiAiTi.1 <- sum(Ri[Ai] * (Ti[Ai] - Oi))/Np

  MiAi.1 <- sum(Ri[Ai] * (Ti[Ai] - Oi))
  MiAi.2 <- sum((1-Ri[Ai]) * (1-Di[Ai]) * (Ti[Ai] - Oi))

  (EMiAiTi - EMiAiOi) == (MiAi.1 + MiAi.2) / Np

  # EM.1 <- beta.est + (OMEGA - beta.est * 2) * INCIDENCE * (1-p) / p
  # EMA.1 <- EM.1 * pA
  # EMAT.1 <- EMA.1 * (muTA - omegaTA)

  EAiRi0Di0Ti <- mean(Ai * (1-Ri) * (1-Di) * Ti)
  EAR0D0T <- INCIDENCE * (1-p) / p * pA * (muTA**2 + varTA - omegaTAstar)

  EAiRi0Di0 <- mean(Ai * (1-Ri) * (1-Di))
  EAR0D0 <- INCIDENCE * (1-p) / p * pA * (muTA - omegaTA)

  EAiRi0Di0Oi <- sum((1-Ri[Ai]) * (1-Di[Ai]) * Oi)/Np
  EAR0D0O <- INCIDENCE * (1-p) / p * pA * (omegaTAstar - (omegaTA**2 + omegaTAvar))

  # browser()
  w2 <- Np

  cov12.1 <- N * p * (1-p) * pr_star
  cov12.2 <- N * p * (1-p) * (1-pB) * beta.est
  cov12 <- N * p * (1-p) * (pr_star - (1-pB) * beta.est)

  return(c(w41=w41, w42=w42, w1=w1, w4=w4, w5=w5, cov45=cov45, var4=var4,
           var41=var41, var42=var42, var43=var43, cov14=cov14,
           cov14.1=cov14.1, cov14.2=cov14.2, cov14.3=cov14.3, cov14.4=cov14.4,
           cov14.intermediate=cov14.intermediate,
           cov14.intermediate.2=cov14.intermediate.2,
           w11=w11, w12=w12,
           EMi.1=EMi.1, # EM.1=EM.1,
           EMiAi.1=EMiAi.1, # EMA.1=EMA.1,
           EMiAiTi.1=EMiAiTi.1, EMAT.1=EMAT.1,
           EAiRi0Di0Ti=EAiRi0Di0Ti, EAR0D0T=EAR0D0T,
           EAiRi0Di0=EAiRi0Di0, EAR0D0=EAR0D0,
           EAiRi0Di0Oi=EAiRi0Di0Oi, EAR0D0O=EAR0D0O,
           w2=w2, cov12.1=cov12.1, cov12.2=cov12.2, cov12=cov12,
           EMiAiTi=EMiAiTi, EMiAiOi=EMiAiOi))
}

sims <- pbmclapply(1:NSIM, one.study, mc.cores=8L)
simdf <- rbindlist(lapply(sims, function(x) as.data.table(t(x))))

cov(simdf$w4, simdf$w5)
mean(simdf$cov45)

var(simdf$w4)
mean(simdf$var4)

var(simdf$w41)
mean(simdf$var41)

var(simdf$w42)
mean(simdf$var42)

cov(simdf$w41, simdf$w42)
mean(simdf$var43)

var(simdf$w41) + var(simdf$w42) - 2 * cov(simdf$w41, simdf$w42)
mean(simdf$var41) + mean(simdf$var42) - 2 * mean(simdf$var43)

# Covariance of W1, W4

cov(simdf$w1, simdf$w4)
mean(simdf$cov14)

cov(simdf$w11, simdf$w41)
mean(simdf$cov14.1)

cov(simdf$w11, simdf$w42)
mean(simdf$cov14.2)

cov(simdf$w12, simdf$w41)
mean(simdf$cov14.3)

cov(simdf$w12, simdf$w42)
mean(simdf$cov14.4)

cov(simdf$w11, simdf$w41) - cov(simdf$w11, simdf$w42) - cov(simdf$w12, simdf$w41) + cov(simdf$w12, simdf$w42)

mean(simdf$cov14.intermediate)
mean(simdf$cov14.intermediate.2)

# Really similar

mean(simdf$EM.1)
mean(simdf$EMi.1)

mean(simdf$EMA.1)
mean(simdf$EMiAi.1)

mean(simdf$EMAT.1)
mean(simdf$EMiAiTi.1)

# Next

mean(simdf$EAiRi0Di0Ti)
mean(simdf$EAR0D0T)

mean(simdf$EAiRi0Di0)
mean(simdf$EAR0D0)

mean(simdf$EAiRi0Di0Oi)
mean(simdf$EAR0D0O)

cov(simdf$w1, simdf$w2)
mean(simdf$cov12)

cov(simdf$w11, simdf$w2)
mean(simdf$cov12.1)

cov(simdf$w12, simdf$w2)
mean(simdf$cov12.2)

