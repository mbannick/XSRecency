set.seed(21)

# PARAMETERS --------
bigT <- 1
tau <- 10
lambda <- 0.03
prevalence <- 0.2
n_sims <- 2
N <- 10
a <- 0.5
b <- 1.0
q <- 0.4
rho <- 0

truephifunc_MDRI <- function(t){
  ###### MDRI at 1 year 130/365, FRR 1.6%: 1.6% + Gamma distribution with mean 130/365-1.6%
  alpha = 3; x = (130/365-0.016)/(1-0.016); beta=alpha/x
  return((1-pgamma(t, shape = alpha, rate = beta))*(1-0.016)+0.016)
}

dat <- generate.raw.data(n_sims=n_sims, n=N, prevalence=prevalence)

sim <- simulate.recent(sim_data=dat, infection.function=infections.con,
                       baseline_incidence=lambda, prevalence=prevalence, rho=0,
                       phi.func=truephifunc_MDRI,
                       summarize=TRUE,
                       ptest.dist=function(n, u) runif(n, a, b),
                       ptest.prob=q, bigT=bigT)
