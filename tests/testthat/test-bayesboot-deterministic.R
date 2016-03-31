context("Deterministic Bayesian bootstrap tests")
library(doParallel)
set.seed(123)

test_that("rudirichlet produces a valid output", {
  rand_mat <- rudirichlet(10, 15)
  expect_true(all(rand_mat >= 0 & rand_mat <= 1))
  expect_equivalent(rowSums(rand_mat), rep(1, 10))
})

# TODO: Why does this pass when using test(), but not when
# checking the package?
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

  # A "stranger" bootstrap analysis with the data being chars. in a list.
  # And the statistc being the most common answer.
  data_list <- list("Yes", "Yes", "No", "Yes", "No", "Yes", "Maybe")
  boot_stat <- function(d) {
    t <- table(as.character(d))
    c(most_common_answer = names(t)[ which.max(t)])
  }
  b5 <- bayesboot(data_list, boot_stat, R = 50, R2 = 20)
  expect_equal(class(b5), c("bayesboot", "data.frame"))
  expect_equal(nrow(b5), 50)
  expect_equal(ncol(b5), 1)

  # Another strange bootstrap with a statistic that outputs NAs
  d <- data.frame(x = 1:15, y = c(1, 2, 3, 4, NA))
  expect_warning({
    b6 <- bayesboot(d, use.weights = TRUE, R = 100, function(d, w) {
      c(weighted.mean(d$x, w), weighted.mean(d$y, w))
    })
  })

  expect_output(print(summary(b1)), ".")
  expect_output(print(summary(b2)), ".")
  expect_output(print(summary(b3)), ".")
  expect_output(print(summary(b4)), ".")
  expect_warning(summary(b5))
  expect_true({
    plot(b1)
    plot(b2)
    plot(b3)
    plot(b4)
    TRUE
  })
  expect_warning(plot(b5))
})

test_that("bayesboot can do paralell processing", {
  library(doParallel)
  library(foreach)
  x <- rnorm(10)
  registerDoParallel(cores = 2)
  b1 <- bayesboot(x, mean, R = 1000, R2 = 1000, .parallel = TRUE)
  expect_equal(class(b1), c("bayesboot", "data.frame"))
  expect_equal(nrow(b1), 1000)
  expect_equal(ncol(b1), 1)
  stopImplicitCluster()
  registerDoParallel(cores = 1)
  stopImplicitCluster()

})

