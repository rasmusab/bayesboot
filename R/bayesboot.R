# Importing some of the base packages which are used in the functions below.
#' @import stats
#' @import grDevices
#' @import graphics
#' @import utils
NULL

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
#' # dataset of size three.
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


#' The Bayesian bootstrap
#'
#' Performs a Bayesian bootstrap and returns a \code{data.frame} with a sample
#' of size \code{R} representing the posterior distribution of the (possibly
#' multivariate) summary \code{statistic}.
#'
#' The summary statistic is a function of the data that represents a feature of
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
#' larger range of statistics (the \code{\link{median}}, for example)
#'
#' For more information regarding this implementation of the Bayesian bootstrap
#' see the blog post
#' \href{http://www.sumsar.net/blog/2015/07/easy-bayesian-bootstrap-in-r/}{Easy
#' Bayesian Bootstrap in R}. For more information about the model behind the
#' Bayesian bootstrap see the blog post
#' \href{http://www.sumsar.net/blog/2015/04/the-non-parametric-bootstrap-as-a-bayesian-model/}{The
#' Non-parametric Bootstrap as a Bayesian Model} and, of course,
#' \href{http://projecteuclid.org/euclid.aos/1176345338}{the original Bayesian
#' bootstrap paper by Rubin (1981)}.
#'
#' @note \itemize{
#' \item  While \code{R} and \code{R2} are set to \code{4000} by
#' default, that should not be taken to indicate that a sample of size 4000 is
#' sufficient nor recommended.
#'
#' \item When using \code{use.weights = FALSE} it is important to use a summary
#' statistic that does not depend on the sample size. That is, doubling the size
#' of a dataset by cloning data should result in the same statistic as when
#' using the original dataset. An example of a statistic that depends on the
#' sample size is the sample standard deviation (that is, \code{\link{sd}}), and
#' when using \code{bayesboot} it would make more sense to use the population
#' standard deviation (as in the example below). }
#'
#' @param data Either a vector or a list, or a matrix or a data.frame with one
#'   datapoint per row. The format of \code{data} should be compatible with the
#'   first argument of \code{statistic}
#' @param statistic A function implementing the summary statistic of interest
#'   where the first argument should take the data. If \code{use.weights = TRUE}
#'   then the second argument should take a vector of weights.
#' @param R The size of the posterior sample from the Bayesian bootstrap.
#' @param R2 When \code{use.weights = FALSE} this is the size of the resample of
#'   the data used to approximate the weighted statistic.
#' @param use.weights When \code{TRUE} the data will be reweighted, like in the
#'   original Bayesian bootstrap. When \code{FALSE} (the default) the
#'   reweighting will be approximated by resampling the data.
#' @param .progress The type of progress bar ("none", "text", "tk", and "win").
#'   See the \code{.progress} argument to \code{\link[plyr]{adply}} in the plyr
#'   package.
#' @param .parallel If \code{TRUE} enables parallel processing. See the
#'   \code{.parallel} argument to \code{\link[plyr]{adply}} in the plyr package.
#' @param ... Other arguments passed on to \code{statistic}
#'
#' @return A \code{data.frame} with \code{R} rows, each row being a draw from
#'   the posterior distribution of the Bayesian bootstrap. The number of columns
#'   is decided by the length of the output from \code{statistic}. If
#'   \code{statistic} does not return a vector or data frame with named values
#'   then the columns will be given the names \code{V1}, \code{V2}, \code{V3},
#'   etc. While the output is a \code{data.frame} it has subclass
#'   \code{bayesboot} which enables specialized \code{\link{summary}} and
#'   \code{\link{plot}} functions for the result of a \code{bayesboot} call.
#'
#' @examples
#'
#' ### A Bayesian bootstrap analysis of a mean ###
#'
#' # Heights of the last ten American presidents in cm (Kennedy to Obama).
#' heights <- c(183, 192, 182, 183, 177, 185, 188, 188, 182, 185);
#' b1 <- bayesboot(heights, mean)
#' # But it's more efficient to use the a weighted statistic.
#' b2 <- bayesboot(heights, weighted.mean, use.weights = TRUE)
#'
#' # The result of bayesboot can be plotted and summarized
#' plot(b2)
#' summary(b2)
#'
#' # It can also be easily post processed.
#' # Here the probability that the mean is > 182 cm.
#' mean( b2[,1] > 182)
#'
#' ### A Bayesian bootstrap analysis of a SD ###
#'
#' # When use.weights = FALSE it is important that the summary statistics
#' # does not change as a function of sample size. This is the case with
#' # the sample standard deviation, so here we have to implement a
#' # function calculating the population standard deviation.
#' pop.sd <- function(x) {
#'   n <- length(x)
#'   sd(x) * sqrt( (n - 1) / n)
#' }
#'
#' b3 <- bayesboot(heights, pop.sd)
#' summary(b3)
#'
#' ### A Bayesian bootstrap analysis of a correlation coefficient ###
#'
#' # Data comparing two methods of measuring blood flow.
#' # From Table 1 in Miller (1974) and used in an example
#' # by Rubin (1981, p. 132).
#' blood.flow <- data.frame(
#'   dye = c(1.15, 1.7, 1.42, 1.38, 2.80, 4.7, 4.8, 1.41, 3.9),
#'   efp = c(1.38, 1.72, 1.59, 1.47, 1.66, 3.45, 3.87, 1.31, 3.75))
#'
#' # Using the weighted correlation (corr) from the boot package.
#' library(boot)
#' b4 <- bayesboot(blood.flow, corr, R = 1000, use.weights = TRUE)
#' hist(b4[,1])
#'
#' ### A Bayesian bootstrap analysis of lm coefficients ###
#'
#' # A custom function that returns the coefficients of
#' # a weighted linear regression on the blood.flow data
#' lm.coefs <- function(d, w) {
#'   coef( lm(efp ~ dye, data = d, weights = w) )
#' }
#'
#' b5 <- bayesboot(blood.flow, lm.coefs, R = 1000, use.weights = TRUE)
#'
#' # Plotting the marginal posteriors
#' plot(b5)
#'
#' # Plotting a scatter of regression lines from the posterior
#' plot(blood.flow)
#' for(i in sample(nrow(b5), size = 20)) {
#'   abline(coef = b5[i, ], col = "grey")
#' }
#' @references
#' Miller, R. G. (1974) The jackknife - a review. \emph{Biometrika},
#' \bold{61(1)}, 1--15.
#'
#' Rubin, D. B. (1981). The Bayesian bootstrap. \emph{The annals of statistics},
#' \bold{9(1)}, 130--134.
#' @export
bayesboot <- function(data, statistic, R = 4000, R2 = 4000, use.weights = FALSE,
                      .progress = "none", .parallel = FALSE, ...) {
  call <- match.call()
  # Pick out the first part of statistic matching a legal variable name,
  # just to be used as a label when plotting later.
  statistic.label <- deparse(substitute(statistic))
  match <- regexpr("^(\\w|\\.|_)*", statistic.label[1])
  statistic.label <- regmatches(statistic.label[1], match)

  # Doing some error checks
  if(length(R) != 1 || !is.wholenumber(R) || R < 1) {
    stop("R should be a single whole number >= 1.")
  }

  if ((NROW(data) == 0) || is.null(NROW(data))) {
    stop("no data in call to 'bayesboot'")
  }

  if(use.weights) {
    if (length(formals(statistic)) < 2) {
      stop("If use.weights == TRUE then statistic should take a weight vector as the second argument.")
    }
    w <- rep(1 / NROW(data), NROW(data))
    tryCatch(stat.result <- statistic(data, w, ...),
      error = function(e) {
        message("Applying the statistic on the original data and with uniform weights resulted in an error")
        stop(e)
      }
    )
    # TODO: Should I add some more checks to stat.result? Like, that it contains no NA, values?
    # or should I maybe do these tests to the final posterior samples and issue a varning if
    # there are NAs, NULLs and similar?

  } else { # use.weights == FALSE
    if(length(R2) != 1 || is.na(R2) || !is.wholenumber(R2) || R2 < 1) {
      stop("If use.weights == FALSE then R2 should be a single whole number >= 1.")
    }
    tryCatch(stat.result <- statistic(data, ...),
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
        statistic(data, weights, ...)
      }
    )

  } else {
    if(is.null(dim(data)) || length(dim(data)) < 2) { # data is a list type of object
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
  attr(boot.sample, "call") <- call
  # Warn if boot.sample contains "non-values".
  col.should.warn <- plyr::laply(boot.sample, function(boot.col) {
    any(is.na(boot.col) |is.nan(boot.col) | is.null(boot.col))
  })
  if(any(col.should.warn)) {
    warning(paste(
      "The sample from bayesboot contains either NAs, NaNs or NULLs.",
      "Make sure that your statistic function only return actual values."))
  }
  boot.sample
}



#' Summarize the result of \code{bayesboot}
#'
#' Summarizes the result of a call to \code{bayesboot} by calculating means, SDs,
#' highest density intervals and quantiles of the posterior marginals.
#'
#' @param object The bayesboot object to summarize.
#' @param cred.mass The probability mass to include in the highest density intervals.
#' @param ... Not used.
#'
#' @return A data frame with three columns: (1) \bold{statistic} the name of the
#'   statistic, (2) \bold{measure} the name of the summarizing measure, and (3)
#'   \bold{value} the value of the summarizing measure.
#'
#' @seealso \code{\link[HDInterval]{hdi}} in the HDInterval package for directly calculating
#'   highest density intervals from \code{bayesboot} objects.
#'
#' @export
summary.bayesboot <- function(object, cred.mass = 0.95, ...) {
  bootsum <- plyr::ldply(seq_len(ncol(object)), function(i) {
    s <- object[,i]
    if(!is.numeric(s)) {
      warning(paste("The statistic", names(object)[i] , "was skipped as",
                    "summary.bayesboot can't handle non-numeric statistics."))
      return(data.frame())
    }
    data.frame(statistic   = names(object)[i],
               measure   = c("mean", "sd", "hdi.low", "hdi.high","q2.5%", "q25%", "median" ,"q75%", "q97.5%"),
               value   = c(mean(s), sd(s), HDInterval::hdi(s, cred.mass), quantile(s, c(0.025, 0.25, 0.5, 0.75, 0.975))))
  })
  attr(bootsum, "statistic.label") <- attr(object, "statistic.label")
  attr(bootsum, "call") <- attr(object, "call")
  attr(bootsum, "R") <- nrow(object)
  attr(bootsum, "cred.mass") <- cred.mass
  class(bootsum) <- c("summary.bayesboot", class(bootsum))
  bootsum
}

#' Print the first number of draws from the Bayesian bootstrap
#'
#' @param x The bayesboot object to print.
#' @param n The number of draws to print.
#' @param ... Not used.
#' @export
print.bayesboot <- function(x, n = 10, ...) {
  cat(paste0("The first ", n," draws (out of ", nrow(x) ,") from the Bayesian bootstrap:\n"))
  cat("\n")
  print(as.data.frame(head(x, n)))
  cat(".. ...\n")
  cat("\n")
  cat("Use summary() to produce a summary of the posterior distribution.\n")
  invisible(x)
}

#' @method print summary.bayesboot
#' @export
print.summary.bayesboot <- function(x, ...) {
  stat.table <- plyr::ddply(x, "statistic", function(s) {
    stats <- s$value
    names(stats) <- s$measure
    stats
  })
  cat("Bayesian bootstrap\n")
  cat("\n")
  cat("Number of posterior draws:", attr(x, "R") , "\n")
  cat("\n")
  if(nrow(x) > 0) {
    hdi.percentage <- paste0(round(attr(x, "cred.mass") * 100), "%")
    cat("Summary of the posterior (with", hdi.percentage,"Highest Density Intervals):\n")
    print(stat.table[,c("statistic","mean", "sd", "hdi.low", "hdi.high")], row.names = FALSE)
    cat("\n")
    cat("Quantiles:\n")
    print(stat.table[,c("statistic", "q2.5%", "q25%", "median" ,"q75%", "q97.5%")], row.names = FALSE)
    cat("\n")
  }
  cat("Call:\n", format(attr(x, "call")))
  invisible(x)
}

#' Coerce to a \code{bayesboot} object
#'
#' This converts an object into a data frame and adds the class
#' \code{bayesboot}. Doing this is only useful in the case you would want to use
#' the \code{plot} and \code{summary} methods for \code{bayesboot} objects.
#'
#' @param object Any object that can be converted to a data frame.
#'
#' @return A \code{data.frame} with subclass \code{bayesboot}.
#' @export
as.bayesboot <- function(object) {
  object <- as.data.frame(object)
  class(object) <- c("bayesboot", class(object))
  if(is.null(attr(object, "statistic.label"))) {
    attr(object, "statistic.label") <- ""
  }
  if(is.null(attr(object, "call"))) {
    attr(object, "call") <- ""
  }
  object
}

#' Plot the result of \code{bayesboot}
#'
#' Produces histograms showing the marginal posterior distributions from a
#' \code{bayesboot} call. Uses the \code{\link{plotPost}} function to produce
#' the individual histograms.
#'
#' @param x The bayesboot object to plot.
#' @param cred.mass the probability mass to include in credible intervals, or
#'   NULL to suppress plotting of credible intervals.
#' @param plots.per.page The maximum numbers of plots per page.
#' @param cex,cex.lab,... Further parameters passed on to
#'   \code{\link{plotPost}}.
#'
#' @export
plot.bayesboot <- function(x, cred.mass = 0.95, plots.per.page = 3, cex = 1.2, cex.lab=1.3, ...) {
  old.devAskNewPage <- devAskNewPage()
  old.par <- par(mfrow = c(min(plots.per.page, ncol(x)) , 1) , mar = c(4.1, 4.1, 0.5, 4.1), mgp = c(2.5, 1, 0))
  on.exit({ # revert graphical parameters
    par(old.par)
    devAskNewPage(old.devAskNewPage)
  })
  n.plots <- 0
  for(i in seq_len(ncol(x))) {
    if(!is.numeric(x[, i])) {
      warning(paste("The statistic", names(x)[i] , "was skipped as",
                    "plot.bayesboot can't handle non-numeric statistics."))
      next
    }
    n.plots <- n.plots + 1
    if(n.plots > plots.per.page) {
      devAskNewPage(TRUE)
    }
    if(ncol(x) == 1 && names(x)[i] == "V1" && attr(x, "statistic.label") != "") {
      # There is only one statistic and it has an uninformative default name
      # so use the begining of the function call instead as a statistic, unless
      # it is empty.
      statistic_name <- attr(x, "statistic.label")
    } else { # use the column name
      statistic_name <- names(x)[i]
    }
    plotPost(x[, i], credMass = cred.mass, xlab = statistic_name, cex = cex, cex.lab = cex.lab, ...)
  }
  invisible(NULL)
}


