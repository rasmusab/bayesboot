context("Bayesian bootstrap")

test_that("rudirichlet produces a valid output", {
  rand_mat <- rudirichlet(10, 15)
  expect_equivalent(rowSums(rand_mat), rep(1, 10))
})
