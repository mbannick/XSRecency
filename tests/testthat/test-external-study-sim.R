context("External study simulation")

test_that("External study simulation", {
  set.seed(10)
  func <- function(t) 1 - pgamma(t, 1, 2)
  expect_warning(phi <- assay.properties.sim(func, bigT=2, tau=12))
  expect_equal(phi$mu_est, 0.4831208, tolerance=1e-6)
  expect_equal(phi$mu_var, 0.0003333987, tolerance=1e-6)
  expect_equal(phi$omega_est, 0.4787885, tolerance=1e-6)
  expect_equal(phi$omega_var, 0.0003097912, tolerance=1e-6)
  expect_equal(phi$beta_est, 0.0007204611, tolerance=1e-6)
  expect_equal(phi$beta_var, 5.186902e-07, tolerance=1e-6)
})

