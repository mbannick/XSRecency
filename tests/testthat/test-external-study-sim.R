context("External study simulation")

test_that("External study simulation", {
  set.seed(10)
  phi <- assay.properties.sim(phi.character.1)
  expect_equal(phi$mu$est, 0.4533345, tolerance=1e-7)
  expect_equal(phi$omega$est, 0.3590966, tolerance=1e-7)
  expect_equal(phi$beta$est, 0.01524064, tolerance=1e-7)
})
