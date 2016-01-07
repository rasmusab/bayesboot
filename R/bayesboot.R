
# Performs a Bayesian bootstrap and returns a sample of size n1 representing the
# posterior distribution of the statistic. Returns a vector if the statistic is
# one-dimensional (like for mean(...)) or a data.frame if the statistic is
# multi-dimensional (like for the coefs. of lm).
# Parameters
#   data      The data as either a vector, matrix or data.frame.
#   statistic A function that accepts data as its first argument and possibly
#             the weights as its second, if use_weights is TRUE.
#             Should return a numeric vector.
#   n1        The size of the bootstrap sample.
#   n2        The sample size used to calculate the statistic each bootstrap draw.
#   use_weights  Whether the statistic function accepts a weight argument or
#                should be calculated using resampled data.
#   weight_arg   If the statistic function includes a named argument for the
#                weights this could be specified here.
#   ...       Further arguments passed on to the statistic function.



#' bayesboot
#'
#' @param data
#' @param statistic
#' @param n1
#' @param n2
#' @param use_weights
#' @param weight_arg
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
bayesboot <- function(data, statistic, n1 = 1000, n2 = 1000 , use_weights = FALSE, weight_arg = NULL, ...) {
  # Draw from a uniform Dirichlet dist. with alpha set to rep(1, n_dim).
  # Using the facts that you can transform gamma distributed draws into
  # Dirichlet draws and that rgamma(n, 1) <=> rexp(n, 1)
  dirichlet_weights <- matrix( rexp(NROW(data) * n1, 1) , ncol = NROW(data), byrow = TRUE)
  dirichlet_weights <- dirichlet_weights / rowSums(dirichlet_weights)

  if(use_weights) {
    stat_call <- quote(statistic(data, w, ...))
    names(stat_call)[3] <- weight_arg
    boot_sample <- apply(dirichlet_weights, 1, function(w) {
      eval(stat_call)
    })
  } else {
    if(is.null(dim(data)) || length(dim(data)) < 2) { # data is a list type of object
      boot_sample <- apply(dirichlet_weights, 1, function(w) {
        data_sample <- sample(data, size = n2, replace = TRUE, prob = w)
        statistic(data_sample, ...)
      })
    } else { # data is a table type of object
      boot_sample <- apply(dirichlet_weights, 1, function(w) {
        index_sample <- sample(nrow(data), size = n2, replace = TRUE, prob = w)
        statistic(data[index_sample, ,drop = FALSE], ...)
      })
    }
  }
  if(is.null(dim(boot_sample)) || length(dim(boot_sample)) < 2) {
    # If the bootstrap sample is just a simple vector return it.
    boot_sample
  } else {
    # Otherwise it is a matrix. Since apply returns one row per statistic
    # let's transpose it and return it as a data frame.
    as.data.frame(t(boot_sample))
  }
}
