context("Data generator")

test_that("Data simulation", {
  data <- generate.raw.data(n=10, n_sims=10, prevalence=0.2)
  expect_equal(length(data$n), 10)
  expect_equal(length(data$n_p), 10)
  expect_equal(length(data$n_n), 10)
})
