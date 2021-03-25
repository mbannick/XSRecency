context("Data generator")

PHI <- function(t) 1 - pgamma(t, 1, 2)

test_that("Data simulator", {
  set.seed(10)
  data <- generate.raw.data(n=100, n_sims=3, prevalence=0.4)
  expect_equal(data$n[1, ], c(100, 100, 100))
  expect_equal(data$n_p[1, ], c(36, 40, 43))
  expect_equal(data$n_n[1, ], c(64, 60, 57))
  data <- simulate.recent(sim_data=data, infection.function=infections.con,
                          phi.func=PHI, baseline_incidence=0.03,
                          prevalence=0.4, rho=NA)
  expect_equal(data$n_r, c(0, 0, 2))
})

test_that("Data simulator longitudinal", {
  set.seed(10)
  data <- generate.raw.data(n=100, n_sims=3, prevalence=0.4, times=c(0, 1, 2))
  expect_equal(data$n[1, ], c(100, 100, 100))
  expect_equal(data$n_p[1, ], c(36, 40, 43))
  expect_equal(data$n_n[1, ], c(64, 60, 57))
  data <- simulate.recent(sim_data=data, infection.function=infections.con,
                          phi.func=PHI, baseline_incidence=0.03,
                          prevalence=0.4, rho=NA)
  expect_equal(data$n[1, ], c(100, 100, 100))
  expect_equal(data$n_p[1, ], c(36, 40, 43))
  expect_equal(data$n_n[1, ], c(64, 60, 57))
  # This won't be the same as the previous test because
  # the seed gets messed up by the runif call in simulate.recent
  expect_equal(data$n_r[1, ], c(1, 0, 2))
})

test_that("Data simulation", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  expect_equal(length(data$n), 10)
  expect_equal(length(data$n_p), 10)
  expect_equal(length(data$n_n), 10)
})

test_that("Basic prevalence", {
  data <- generate.raw.data(n=10000, n_sims=10000, prevalence=0.2)
  expect_equal(0.2, mean(data$n_p/data$n), tolerance=1e-4)
})

test_that("Simulate recent, constant", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  data <- simulate.recent(sim_data=data, infection.function=infections.con,
                          phi.func=PHI, baseline_incidence=0.03,
                          prevalence=0.2, rho=NA)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate recent, linear", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  data <- simulate.recent(sim_data=data, infection.function=infections.lin,
                          phi.func=PHI, baseline_incidence=0.03,
                          prevalence=0.2, rho=0.1)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate recent, exponential", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  data <- simulate.recent(sim_data=data, infection.function=infections.exp,
                          phi.func=PHI, baseline_incidence=0.03,
                          prevalence=0.2, rho=0.01)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate full data", {
  data <- generate.data(n=10, n_sims=10,
                        infection.function=infections.exp,
                        phi.func=PHI, baseline_incidence=0.03,
                        prevalence=0.2, rho=0.01)
  expect_equal(length(data$n_r), 10)
})

test_that("Simulate full data unit record", {
  set.seed(10)
  data <- generate.data(n=10, n_sims=10,
                        infection.function=infections.con,
                        phi.func=PHI, baseline_incidence=0.03,
                        prevalence=0.2, rho=NA, summarize=F, times=c(0, 1))
  expect_equal(nrow(data), 10*10*2)
  expect_equal(length(unique(data$rpos)), 3)
  expect_equal(length(unique(data$pos)), 2)
  expect_equal(nrow(data[time == 1]), 10*10)
  expect_equal(nrow(data[time == 0]), 10*10)
  expect_equal(any(is.na(data[pos == 1, itime])), FALSE)
  expect_equal(any(is.na(data[pos == 1, rpos])), FALSE)
  expect_equal(any(is.na(data[pos == 1, itime])), FALSE)
  neg.num <- nrow(data[pos == 0])
  expect_equal(nrow(data[is.na(itime)]), neg.num)
  expect_equal(nrow(data[is.na(rpos)]), neg.num)
  expect_equal(nrow(data[is.na(probs)]), neg.num)
})

