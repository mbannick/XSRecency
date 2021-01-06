context("Estimators")

test_that("Snapshot estimator", {
  set.seed(10)
  data <- generate.data(n=10000, n_sims=100,
                        infection.function=c.infections,
                        phi.func=phi.character.1, baseline_incidence=0.05,
                        prevalence=0.3, rho=NA, frr=0, mdri=142)
  snap <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                       n=data$n, mu_sim=list(est=142/365.25, var=0))
  expect_equal(round(mean(snap$est), 2), 0.05)
})

test_that("Snapshot estimator plus sim mu", {
  set.seed(10)
  assay <- assay.properties.sim(phi.character.1, frr=0)
  data <- generate.data(n=10000, n_sims=100,
                        infection.function=c.infections,
                        phi.func=phi.character.1, baseline_incidence=0.05,
                        prevalence=0.3, rho=NA, frr=0, mdri=142)
  snap <- get.snapshot(n_r=data$n_r, n_n=data$n_n, n_p=data$n_p,
                       n=data$n, mu_sim=list(est=142/365.25, var=0))
  expect_equal(round(mean(snap$est), 2), 0.05)
})
