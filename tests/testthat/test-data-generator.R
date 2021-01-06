context("Data generator")

test_that("Data simulation", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  expect_equal(length(data$n), 10)
  expect_equal(length(data$n_p), 10)
  expect_equal(length(data$n_n), 10)
})

test_that("Simulate recent, constant", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  data <- simulate.recent(sim_data=data, infection.function=c.infections,
                          phi.func=phi.character.1, baseline_incidence=0.03,
                          prevalence=0.2, rho=NA)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate recent, linear", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  data <- simulate.recent(sim_data=data, infection.function=l.infections,
                          phi.func=phi.character.1, baseline_incidence=0.03,
                          prevalence=0.2, rho=0.1)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate recent, exponential", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  data <- simulate.recent(sim_data=data, infection.function=e.infections,
                          phi.func=phi.character.2, baseline_incidence=0.03,
                          prevalence=0.2, rho=0.01, frr=0)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate full data", {
  data <- generate.data(n=10, n_sims=10,
                        infection.function=e.infections,
                        phi.func=phi.character.2, baseline_incidence=0.03,
                        prevalence=0.2, rho=0.01, frr=0)
  expect_equal(length(data$n_r), 10)
})
