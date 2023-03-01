
rm(list=ls())
library(XSRecency)

INC <- 0.032
PREV <- 0.29

N <- 5000
NP <- PREV*N
NN <- N - NP

OMEGA <- 0.25
OMEGAVAR <- 0.0002
BIGT <- 2
MUTB <- 3
VARTB <- 1/3

PRS <- function(pB, beta) INC * (1-PREV) / PREV * (OMEGA - beta * (BIGT - pB * MUTB)) +
  beta * (1-pB)
PR <- function(beta) INC * (1-PREV) / PREV * (OMEGA - beta * BIGT) + beta

PRS(0.5, 0.02) - (PR(0.02) - 0.02 * 0.5 * (1 - INC*(1-PREV)/PREV*MUTB))

vars <- function(pB, beta){
  betavar <- beta * (1-beta) / 1500
  betavar <- 0

  NR <- NP * PR(beta)
  NRPT <- NP * PRS(pB, beta)

  est1 <- adjusted.estimate(
    n_n=NN,
    n_r=NR,
    n_p=NP,
    omega=OMEGA,
    beta=beta,
    big_T=BIGT
  )
  v1 <- variance(
    n_n=NN,
    n_r=NR,
    n_p=NP,
    n=N,
    omega=OMEGA,
    omega_var=OMEGAVAR,
    beta=beta,
    beta_var=betavar,
    big_T=BIGT
  )
  est2 <- adjusted.estimate.pt(
    n_r_pt=NRPT,
    n_n=NN,
    n_p=NP,
    omega=OMEGA,
    beta=beta,
    big_T=BIGT,
    num_beta=NP * (1-pB),
    den_beta=NP * pB * MUTB,
    den_omega=0
  )
  v2s <- variance.pt(
    n_n=NN,
    n_r_pt=NRPT,
    n_r=NR,
    n_p=NP,
    n=N,
    omega=OMEGA,
    omega_var=OMEGAVAR,
    beta=beta,
    beta_var=betavar,
    big_T=BIGT,
    r_TA=0,
    r_TAprime=0,
    r_TAstar=0,
    omega_TA=0,
    omega_TA_var=0,
    omega_TAstar=0,
    omega_TA2=0,
    mu_TA=0,
    var_TA=0,
    mu_TB=MUTB,
    var_TB=VARTB,
    p_A=0,
    p_B=pB)
  v2 <- v2s$logvars
  r1 <- var.log.to.var(estimate=est1, variance=v1)
  r2 <- var.log.to.var(estimate=est2, variance=v2)
  return(r2/r1)
  # return(c(r1, r2, est1, est2, v1, v2, NR, NRPT, v2s$components_est$VW5))
  # return(c(v1, v$logvars))
}



library(ggplot2)
ALLOWED_BETAS <- OMEGA / BIGT # otherwise negative denominator in estimator
pbs <- seq(0, 1, by=0.01)
betas <- seq(0, 0.014, by=0.001)
grid <- expand.grid(pB=pbs, beta=betas)
vals <- mapply(vars, pB=grid$pB, beta=grid$beta)
grid$rvar <- vals

ggplot(data=grid, aes(x=pB, y=rvar, color=beta, group=beta)) +
  geom_point() +
  geom_line() +
  geom_abline(intercept=1.0, slope=0)

PR(beta=0.014)
PRS(pB=0.5, beta=0.014)
