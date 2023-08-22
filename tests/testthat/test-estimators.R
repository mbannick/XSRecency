context("Estimators")
library(parallel)
library(mcreplicate)
library(data.table)

set.seed(10)

# SET PARAMETERS
PHI1 <- function(t) (1 - pgamma(t, 1, 2.7))*(t <= 2) + 0
PHI2 <- function(t) (1 - pgamma(t, 1, 2.7))
N <- 5000
NSIMS <- 5000

test_that("Snapshot estimator", {

  MU <- integratePhi(PHI1, 12)
  set.seed(10)

  get_one <- function(){
    sim <- simCrossSect(
      phi.func=PHI1,
      incidence_type="constant",
      prevalence=0.29,
      baseline_incidence=0.032)

    data <- sim(N)

    snap <- estSnapshot(
      n_r=sum(data$ri, na.rm=T),
      n_n=sum(!data$di),
      n_p=sum(data$di),
      n=N,
      mu=MU, mu_var=0)

    return(snap$est)
  }
  ests <- mc_replicate(n=NSIMS, expr=get_one(), mc.cores=2)

  expect_equal(round(mean(ests), digits=3), 0.032)
})

test_that("Adjusted estimator", {

  set.seed(10)
  MDRI <- integratePhi(PHI2, 2)

  get_one <- function(){
    sim <- simCrossSect(
      phi.func=PHI2,
      incidence_type="constant",
      prevalence=0.29,
      baseline_incidence=0.032)

    data <- sim(N)

    adj <- estAdjusted(
      n_r=sum(data$ri, na.rm=T),
      n_n=sum(!data$di),
      n_p=sum(data$di),
      n=N,
      omega=MDRI, omega_var=0,
      beta=0, beta_var=0,
      big_T=2)

    return(adj$est)
  }
  ests <- mc_replicate(n=NSIMS, expr=get_one(), mc.cores=detectCores())

  expect_equal(round(mean(ests), digits=3), 0.032)
})

test_that("Enhanced estimator", {

  get_one <- function(){

    phidat <- XSRecency:::assay.nsim.pt(
      n_sims=1, phi.func=PHI2, tau=12, bigT=2,
      ext_FRR=FALSE, ext_df=NULL,
      max_FRR=NULL
    )
    phidat <- phidat$studies
    colnames(phidat) <- c("id", "ui", "ri")
    phidat <- phidat[phidat$ui > 0,]

    set.seed(10)
    sim <- simCrossSect(
      phi.func=PHI2,
      incidence_type="constant",
      prevalence=0.29,
      baseline_incidence=0.032)

    df <- sim(N*2)

    sim.pt <- simPriorTests(
      ptest.dist=function(u) runif(1, 0, 4),
      ptest.prob=function(u) 0.5)

    ptdf <- sim.pt(df[df$di == 1,])

    estimate <- estEnhanced(
      n_p=nrow(ptdf),
      n=nrow(df),
      ptdf=ptdf,
      beta=0,
      beta_var=0,
      big_T=2,
      phidat=phidat,
      use_geese=TRUE,
      formula="ri ~ poly(log(ui), degree=2, raw=T)",
      family=binomial(link="logit"),
      min_dt=TRUE
    )

    return(estimate)
  }
  get_one()
  ests <- mc_replicate(n=50, expr=get_one(), mc.cores=detectCores())
  ests <- unlist(ests["est",])
  expect_equal(round(mean(ests), digits=3), 0.032)

})
