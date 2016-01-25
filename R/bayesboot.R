#' Produce random draws from a uniform Dirichlet distribution
#'
#' \code{rudirichlet} produces \code{n} draws from a \code{d}-dimensional
#' uniform Dirichlet distribution. Here "uniform" implies that any combination
#' of values on the support of the distribution is equally likely, that is, the
#' \eqn{\alpha} parameters to the Dirichlet distribution are all set to 1.0.
#'
#' In the context of the Bayesian bootstrap \code{rudirichlet} is used to
#' produces the bootstrap weights. Therefore, \code{rudirichlet} can be used if
#' you directly want to generate Bayesian bootstrap weights.
#'
#'
#' @param n the number of draws.
#' @param d the dimension of the Dirichlet distribution.
#'
#' @return An \code{n} by \code{d} matrix.
#' @export
#'
#' @examples
#' set.seed(123)
#' rudirichlet(2, 3)
#' # Should produce the following matrix:
#' #       [,1]   [,2]   [,3]
#' # [1,] 0.30681 0.2097 0.4834
#' # [2,] 0.07811 0.1390 0.7829
#'
#' # The above could be seen as a sample of two Bayesian bootstrap weights for a
#' dataset of size three.
rudirichlet <- function(n, d) {
  # Using the facts that you can transform gamma distributed draws into
  # Dirichlet draws and that rgamma(n, 1) <==> rexp(n, 1)
  # See here for explanation:
  # https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation#Gamma_distribution
  rexp.mat <- matrix( rexp(d * n, 1) , nrow = n, ncol = d)
  dirichlet.weights <- rexp.mat / rowSums(rexp.mat)
  dirichlet.weights
}

# From the examples on the help page for base::is.integer
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - round(x)) < tol
}

# Performs a Bayesian bootstrap and returns a sample of size R representing the
# posterior distribution of the statistic. Returns a vector if the statistic is
# one-dimensional (like for mean(...)) or a data.frame if the statistic is
# multi-dimensional (like for the coefs. of lm).
# Parameters
#   data      The data as either a vector, matrix or data.frame.
#   statistic A function that accepts data as its first argument and possibly
#             the weights as its second, if use_weights is TRUE.
#             Should return a numeric vector.
#   R        The size of the bootstrap sample.
#   R2       The sample size used to calculate the statistic each bootstrap draw.
#   use.weights  Whether the statistic function accepts a weight argument or
#                should be calculated using resampled data.
#   ...       Further arguments passed on to the statistic function.



#' Bayesian bootstrap
#'
#' Performs a Bayesian bootstrap and returns a \code{data.frame} with a sample
#' of size \code{R} representing the posterior distribution of the (possibly
#' multivariate) statistic.
#'
#' The statistic is a function of the data that represents a feature of
#' interest, where a typical statistic is the mean. In \code{bayesboot} it is
#' most efficient to define the statistic as a function taking the data as the
#' first argument and a vector of weights as the second argument. An example of
#' such a function is \code{\link{weighted.mean}}. Indicate that you are using a
#' statistic defined in this way by setting \code{use.weights = TRUE}.
#'
#' It is also possible to define the statistic as a function only taking data
#' (and no weights) by having \code{use.weights = FALSE} (the default). This
#' will, for each of the \code{R} Bayesian bootstrap draws, give a resampled
#' version of the \code{data} of size \code{R2} to \code{statistic}. This will
#' be much slower than using \code{use.weights = TRUE} but will work with a
#' larger range of statistics (\code{\link{median}}, for example)
#'
#' \emph{Note 1}: While \code{R} and \code{R2} are set to \code{1000} by default,
#' that should not be taken to indicate that a sample of size 1000 is sufficient nor recommended.
#'
#' \emph{Note 2}: When using \code{use.weights = FALSE} it is important to use
#' a statistic that does not depend on the sample size. That is, doubling the size of a dataset by cloning data
#' should result in the same statistic as the original dataset. An example of a statistic that depends on
#' the sample size is the sample standard deviation (that is, \code{\link{sd}}), and when using \code{bayesboot} it
#' would make more sense to use the population standard deviation (see examples below).
#'
#'
#' @param data Either a vector or a list, or a matrix or a data.frame with one
#'   datapoint per row. The format of \code{data} should be compatible with the
#'   first argument of \code{statistic}
#' @param statistic A function implementing the statistic of interest where the
#'   first argument should take the data. If \code{use.weights = TRUE} then the
#'   second argument should take a vector of weights.
#' @param R The size of the posterior sample from the Bayesian bootstrap.
#' @param R2 When \code{use.weights = FALSE} this is the size of the resample of
#'   the data used to approximate the weighted statistic.
#' @param use.weights When \code{TRUE} the data will be reweighted, like in the
#'   original Bayesian bootstrap. When \code{FALSE} (the default) the reweighting
#'   will be approximated by resampling the data.
#' @param .progress The type of progress bar. See the \code{.progress} argument
#'   to \code{\link[plyr]{adply}} in the plyr package.
#' @param .parallel If \code{TRUE} enables parallel processing. See the
#'   \code{.parallel} argument to \code{\link[plyr]{adply}} in the plyr package.
#' @param ... Other arguments passed on to \code{statistic}
#'
#' @return A \code{data.frame} with \code{R} rows, each row being a draw from
#'   the posterior distribution of the Bayesian bootstrap. The number of columns
#'   is decided by the length of the output from \code{statistic}. While the
#'   output is a \code{data.frame} it has subclass \code{bayesboot} which
#'   enables specialized \code{\link{summary}} and \code{\link{plot}} functions
#'   for the result of a \code{bayesboot} call.
#'
#' @examples
#' # A function calculating the population standard deviation.
#' pop.sd <- function(x) {
#'   n <- length(x)
#'   sd(x) * sqrt( (n - 1) / n)
#' }
#' @export
bayesboot <- function(data, statistic, R = 1000, R2 = 1000, use.weights = FALSE,
                      .progress = "none", .parallel = FALSE, ...) {
  # Pick out the first part of statistic matching a legal variable name,
  # just to be used as a label when plotting later.
  statistic.label <- deparse(substitute(statistic))
  match <- regexpr("^(\\w|\\.|_)*", "function <- fun")
  statistic.label <- regmatches(statistic.label, match)

  # Doing some error checks
  if(length(R) != 1 || !is.wholenumber(R) || R < 1) {
    stop("R should be a single whole number >= 1.")
  }

  if ((NROW(data) == 0) || is.null(NROW(data))) {
    stop("no data in call to 'bayesboot'")
  }

  if(use.weights) {
    if(length(R2) != 1 || !is.wholenumber(R2) || R2 < 1) {
      stop("If use.weights == TRUE then R2 should be a single whole number >= 1.")
    }
    if (length(formals(statistic)) < 2) {
      stop("If use.weights == TRUE then statistic should take a weight vector as the second argument.")
    }
    w <- rep(1 / NROW(data), NROW(data))
    tryCatch(stat.result <- statistic(data, w),
      error = function(e) {
        message("Applying the statistic on the original data and with uniform weights resulted in an error")
        stop(e)
      }
    )
  } else { # use.weights == FALSE
    tryCatch(stat.result <- statistic(data),
      error = function(e) {
        message("Applying the statistic on the original data resulted in an error")
        stop(e)
      }
    )
  }

  if(! (is.atomic(stat.result) || is.data.frame(stat.result) || is.matrix(stat.result)) &&
        NROW(stat.result) != 1) {
    stop(paste(
    "Applying the statistic on the original data should return a vector, or a",
     "data.frame or matrix with one row, but an object with the following",
     "structure was returned instead:\n",
    paste(capture.output(str(stat.result)), collapse = "\n")
    ))
  }

  dirichlet.weights <- rudirichlet(R, NROW(data))

  if(use.weights) {
    boot.sample <- plyr::adply(
      dirichlet.weights, 1, .progress = .progress, .parallel = .parallel, .id = NULL,
      .fun = function(weights) {
        statistic(data, weights)
      }
    )

  } else {
    if(is.atomic(x) || is.list(x)) { # data is a list type of object
      boot.sample <- plyr::adply(
        dirichlet.weights, 1, .progress = .progress, .parallel = .parallel, .id = NULL,
        .fun = function(weights) {
          data.sample <- sample(data, size = R2, replace = TRUE, prob = weights)
          statistic(data.sample, ...)
        }
      )
    } else { # assume data can be subsetted like a matrix or data.frame
      boot.sample <- plyr::adply(
        dirichlet.weights, 1, .progress = .progress, .parallel = .parallel, .id = NULL,
        .fun = function(weights) {
          index.sample <- sample(nrow(data), size = R2, replace = TRUE, prob = weights)
          statistic(data[index.sample, ,drop = FALSE], ...)
        }
      )
    }
  }
  class(boot.sample) <- c("bayesboot", class(boot.sample))
  attr(boot.sample, "statistic.label") <- statistic.label
  boot.sample
}


