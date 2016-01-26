# This file contains some tests that are non-deterministic and that might
# occasionally fail, even if everything is OK. But they shouldn't fail too often...

context("Non-deterministic Bayesian bootstrap tests")


test_that("rudirichlet produces a valid output", {
  rand_mat <- rudirichlet(10, 15)
  expect_equivalent(rowSums(rand_mat), rep(1, 10))
})
