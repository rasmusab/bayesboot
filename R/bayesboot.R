#' Produce random draws from a uniform Dirichlet distribution
#'
#' \code{rudirichlet} produces \code{n} draws from a \code{d}-dimensional uniform Dirichlet distribution.
#' Here "uniform" implies that any combination of values on the support of the distribution is equally likely,
#' that is, the \eqn{\alpha} parameters to the Dirichlet distribution are all set to 1.0.
#'
#'
#' @param n the number of draws.
#' @param d the dimension of the Dirichlet distribution.
#'
#' @return An \code{n} by \code{d} matrix.
#' @export
#'
#' @examples
#' rudirichlet(3, 3)
#' # Should produce the following matrix, but with different numbers, of course:
#' #        [,1]     [,2]   [,3]
#' # [1,] 0.7109 0.001012 0.2881
#' # [2,] 0.2112 0.491784 0.2970
rudirichlet <- function(n, d) {
  # Using the facts that you can transform gamma distributed draws into
  # Dirichlet draws and that rgamma(n, 1) <==> rexp(n, 1)
  # See here for explanation:
  # https://en.wikipedia.org/wiki/Dirichlet_distribution#Random_number_generation#Gamma_distribution
  rexp_mat <- matrix( rexp(d * n, 1) , ncol = d, byrow = TRUE)
  dirichlet_weights <- rexp_mat / rowSums(rexp_mat)
  dirichlet_weights
}


# Performs a Bayesian bootstrap and returns a sample of size n1 representing the
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
#   weight.arg   If the statistic function includes a named argument for the
#                weights this could be specified here.
#   ...       Further arguments passed on to the statistic function.
bayesboot <- function(data, statistic, R = 1000, R2 = 1000, use.weights = FALSE, weight.arg = NULL, ...) {
  # Pick out the first part of statistic matching a legal variable name,
  # just to be used as a label when plotting later.
  statistic.label <- deparse(substitute(statistic))
  match <- regexpr("^(\\w|\\.|_)*", "function <- fun")
  statistic.label <- regmatches(statistic.label, match)

  # Todo: Add a lot of error checks...

  dirichlet_weights <- rudirichlet(R, NROW(data))

  if(use.weights) {
    stat_call <- quote(statistic(data, weights, ...))
    names(stat_call)[3] <- weight.arg
    boot_sample <- adply(dirichlet_weights, 1, .id = NULL, function(weights) {
      eval(stat_call)
    })
  } else {
    if(is.atomic(x) || is.list(x)) { # data is a list type of object
      boot_sample <- adply(dirichlet_weights, 1, .id = NULL, function(weights) {
        data_sample <- sample(data, size = R2, replace = TRUE, prob = weights)
        statistic(data_sample, ...)
      })
    } else { # assume data can be subsetted like a matrix or data.frame
      boot_sample <- adply(dirichlet_weights, 1, .id = NULL, function(weights) {
        index_sample <- sample(nrow(data), size = R2, replace = TRUE, prob = weights)
        statistic(data[index_sample, ,drop = FALSE], ...)
      })
    }
  }
  class(boot_sample) <- c("bayesboot", class(boot_sample))
  attr(boot_sample, "statistic.label") <- statistic.label
  boot_sample
}
