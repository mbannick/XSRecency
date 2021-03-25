context("External study simulation")

test_that("External study simulation", {
  set.seed(200)
  func <- function(t) 1 - pgamma(t, 1, 2)
  expect_warning(phi <- assay.properties.nsim(1, func, bigT=2, tau=12,
                                              integrate.FRR=FALSE))
  expect_equal(phi$mu_est, 0.4902688, tolerance=1e-6)
  expect_equal(phi$mu_var, 4.238794e-07, tolerance=1e-6)
  expect_equal(phi$omega_est, 0.4793774, tolerance=1e-6)
  expect_equal(phi$omega_var, 3.745049e-07, tolerance=1e-6)
  expect_equal(phi$beta_est, 0.0007204611, tolerance=1e-6)
  expect_equal(phi$beta_var, 5.186902e-07, tolerance=1e-6)
})

