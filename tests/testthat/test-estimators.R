context("Estimators")

PHI <- function(t) 1 - pgamma(t, 1, 2.7)

test_that("Snapshot estimator", {
  set.seed(10)
  data <- generate.data(n=10000, n_sims=100,
                        infection.function=c.infections,
                        phi.func=PHI, baseline_incidence=0.05,
                        prevalence=0.3, rho=NA)
  snap <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                       n=data$n, mu=142/365.25, mu_var=0)
  expect_equal(round(mean(snap$est), 2), 0.05)
})

test_that("Snapshot estimator plus sim mu", {
  set.seed(10)
  expect_warning(assay <- assay.properties.nsim(100, PHI, bigT=2, tau=12))
  data <- generate.data(n=10000, n_sims=100,
                        infection.function=c.infections,
                        phi.func=PHI, baseline_incidence=0.05,
                        prevalence=0.3, rho=NA)
  snap <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                       n=data$n, mu=assay$mu_est, mu_var=assay$mu_var)
  expect_equal(round(mean(snap$est), 2), 0.05)
})

test_that("Adjusted estimator", {
  set.seed(10)
  data <- generate.data(n=10000, n_sims=100,
                        infection.function=c.infections,
                        phi.func=PHI, baseline_incidence=0.05,
                        prevalence=0.3, rho=NA)
  adj <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                       n=data$n, omega=142/365.25, omega_var=0,
                       beta=0, beta_var=0, big_T=3)
  expect_equal(round(mean(adj$est), 2), 0.05)
})

test_that("Adjusted estimator plus sim omega, beta", {
  set.seed(10)
  expect_warning(assay <- assay.properties.nsim(n_sims=100, PHI, bigT=2, tau=12))
  data <- generate.data(n=10000, n_sims=100,
                        infection.function=c.infections,
                        phi.func=PHI, baseline_incidence=0.05,
                        prevalence=0.3, rho=NA)
  adj <- get.adjusted(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                      n=data$n, omega=assay$omega_est,
                      omega_var=assay$omega_var,
                      beta=assay$beta_est, beta_var=assay$beta_var, big_T=3)
  expect_equal(round(mean(adj$est), 2), 0.05)
})
