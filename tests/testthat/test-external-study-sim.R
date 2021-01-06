context("External study simulation")

test_that("External study simulation", {
  set.seed(10)
  phi <- assay.properties.sim(phi.character.1)
  expect_equal(phi$mu$est, 0.4533345, tolerance=1e-7)
  expect_equal(phi$omega$est, 0.3590966, tolerance=1e-7)
  expect_equal(phi$beta$est, 0.01524064, tolerance=1e-7)
})

test_that("External study simulation 2", {
  set.seed(10)
  phi <- assay.properties.sim(phi.character.2)
  expect_equal(phi$mu$est, 0.8880401, tolerance=1e-7)
  expect_equal(phi$omega$est, 0.3591058, tolerance=1e-7)
  expect_equal(phi$beta$est, 0.01524064, tolerance=1e-7)
})

test_that("External study simulation 3", {
  set.seed(10)
  phi <- assay.properties.sim(phi.character.1, frr=0)
  expect_equal(phi$mu$est, 0.369737, tolerance=1e-7)
  expect_equal(phi$omega$est, 0.349273, tolerance=1e-7)
  expect_equal(phi$beta$est, 0.01524064, tolerance=1e-7)
})
