context("Deterministic Bayesian bootstrap tests")
set.seed(123)

test_that("rudirichlet produces a valid output", {
  rand_mat <- rudirichlet(10, 15)
  expect_true(all(rand_mat >= 0 & rand_mat <= 1))
  expect_equivalent(rowSums(rand_mat), rep(1, 10))
})

test_that("bayesboot produces a valid output", {
  x <- rnorm(10)

  b1 <- bayesboot(x, mean, R = 100, R2 = 90, use.weights = FALSE)
  expect_equal(class(b1), c("bayesboot", "data.frame"))
  expect_equal(nrow(b1), 100)
  expect_equal(ncol(b1), 1)

  b2 <- bayesboot(x, weighted.mean, R = 50, R2 = NULL, use.weights = TRUE)
  expect_equal(class(b2), c("bayesboot", "data.frame"))
  expect_equal(nrow(b2), 50)
  expect_equal(ncol(b2), 1)

  d <- data.frame(x = 1:10, y = rnorm(10))

  boot_stat <- function(d) {
    coef(lm(y ~ x, data = d))
  }
  b3 <- bayesboot(d, boot_stat, R = 75, R2 = 1000, use.weights = FALSE)
  expect_equal(class(b3), c("bayesboot", "data.frame"))
  expect_equal(nrow(b3), 75)
  expect_equal(ncol(b3), 2)

  boot_stat <- function(d, w) {
    coef(lm(y ~ x, data = d, weights = w))
  }
  b4 <- bayesboot(d, boot_stat, R = 130, use.weights = TRUE)
  expect_equal(class(b4), c("bayesboot", "data.frame"))
  expect_equal(nrow(b4), 130)
  expect_equal(ncol(b4), 2)

  expect_true({
    summary(b1)
    summary(b2)
    summary(b3)
    summary(b4)

    #plot(b1)
    #plot(b2)
    #plot(b3)
    #plot(b4)
    TRUE
  })
})
