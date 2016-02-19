# This file contains some tests that are non-deterministic and that might
# occasionally fail, even if everything is OK. But they shouldn't fail too often...

context("Non-deterministic Bayesian bootstrap tests")

x = runif(5)
blood.flow <- data.frame(
  dye = c(1.15, 1.7, 1.42, 1.38, 2.80, 4.7, 4.8, 1.41, 3.9),
  efp = c(1.38, 1.72, 1.59, 1.47, 1.66, 3.45, 3.87, 1.31, 3.75))

test_that("The weight based and the resampling based bayesboot does the same thing", {
  skip_on_cran() #As these tests might occasionally fail.
  b1 <- bayesboot(x, mean, R = 10000, R2 = 1000)
  b2 <- bayesboot(x, weighted.mean, R = 10000)
  q1 <- quantile(b1$V1, c(0.1, 0.5, 0.9))
  q2 <- quantile(b2$V1, c(0.1, 0.5, 0.9))
  # Check that some quantiles are roughly the same
  expect_true(all(abs(q1 - q2) < 0.01))
})

test_that("bayesboot replicates the correlation anlalysis in Rubin (1981)", {
  skip_on_cran() #As these tests might occasionally fail.
  library(boot)
  target_q <- c("0.1" = 0.8962, "0.5" = 0.9519, "0.9" = 0.9788)
  b3 <- bayesboot(blood.flow, boot::corr, R = 10000, use.weights = TRUE)
  b4 <- bayesboot(blood.flow, function(x) { cor(x[,1], x[,2]) }, R = 10000, R2 = 1000)
  q3 <- quantile(b3$V1, c(0.1, 0.5, 0.9))
  q4 <- quantile(b4$V1, c(0.1, 0.5, 0.9))
  expect_true(all(abs(q3 - target_q) < 0.005))
  expect_true(all(abs(q4 - target_q) < 0.005))
})
