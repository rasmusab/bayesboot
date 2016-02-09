#' Highest Density Interval
#'
#' Calculate the highest density interval (HDI) for a given probability mass,
#' see Details. The function is generic, with methods for a range of input
#' objects.
#'
#' The HDI is the interval which contains the required mass such that all
#' points within the interval have a higher probability density than points
#' outside the interval. When applied to a posterior probability density, it is
#' often known as the Highest Posterior Density (HPD).
#'
#' \figure{HDIskew.jpg}
#'
#' In contrast, a symmetric density interval defined by (eg.) the 10\% and 90\%
#' quantiles may include values with lower probability than those excluded.
#'
#' For a unimodal distribution, the HDI is the narrowest interval containing
#' the specified mass, and the \code{hdi} function actually returns the
#' narrowest interval. This does not always work properly for multimodel
#' densities, where the HDI may be discontinuous (the horizontal black line in
#' the Figure below). The single interval returned by \code{hdi} (the blue
#' line) may incorrectly include values between the modes with low probability
#' density. The \code{density} method with \code{allowSplit = TRUE} gives
#' separate limits for discontinuous HDIs.
#'
#' \figure{HDIbimodal.jpg}
#'
#' The default method expects a vector representing draws from the target
#' distribution, such as is produced by an MCMC process. Missing values are
#' silently ignored; if the vector has no non-missing values, NAs are returned.
#'
#' The matrix and data frame methods expect an object with vectors of the above
#' type for each parameter in columns. The result is a matrix with parameters
#' in columns, and rows with the upper and lower limits of the HDI.
#'
#' The mcmc.list method expects an object of type \code{mcmc.list} as defined
#' in package \pkg{coda}.
#'
#' None of the above use interpolation: the values returned correspond to
#' specific values in the data object. Results thus depend on the random draws,
#' and will be unstable if few values are provided. For a 95\% HDI, 10,000
#' independent draws are recommended; a smaller number will be adequate for a
#' 80\% HDI, many more for a 99\% HDI.
#'
#' The function method requires the name for the inverse cumulative density
#' function (ICDF) of the distribution; standard R functions for this have a
#' \code{q-} prefix, eg. \code{qbeta}. Arguments required by the ICDF must be
#' specified by their (abbreviated) names; see the examples.
#'
#' @param object an object specifying the target distribution; see Details.
#' @param credMass a scalar [0, 1] specifying the mass within the credible
#' interval.
#' @param tol the desired accuracy; see \code{\link{optimize}}; default is
#' 1e-8.
#' @param allowSplit if TRUE and the HDI is discontinuous, the beginning and
#' end of each segment are returned; see Value.
#' @param \dots named parameters to be passed to other methods; see Examples.
#' @aliases hdi hdi.default hdi.function hdi.matrix hdi.data.frame hdi.density
#'   hdi.mcmc.list
#' @return a vector of length 2 or a 2-row matrix with the lower and upper
#' limits of the HDI, with an attribute "credMass".
#'
#' The \code{density} method with \code{allowSplit = TRUE} produces a matrix
#' with a row for each component of a discontinuous HDI and columns for begin
#' and end. It has an additional attribute "height" giving the probability
#' density at the limits of the HDI.
#' @note The origin of this function is
#' \href{https://cran.r-project.org/web/packages/BEST/index.html}{the BEST package}.
#'
#' @author Mike Meredith. Code for \code{hdi.function} based on \code{hpd} by
#' Greg Snow, corrected by John Kruschke.
#' @references Kruschke, J. K. (2015) \emph{Doing Bayesian data analysis, second
#' edition: A tutorial with R, JAGS, and Stan.} Waltham, MA: Academic Press /
#' Elsevier.
#'
#' Kruschke, J. K. (2013) Bayesian estimation supersedes the t test.
#' \emph{Journal of Experimental Psychology: General}, \bold{142(2)}, 573.
#' @examples
#'
#' # for a vector:
#' tst <- rgamma(1e5, 2.5, 2)
#' hdi(tst)
#' hdi(tst, credMass=0.8)
#' # For comparison, the symmetrical 80% CrI:
#' quantile(tst, c(0.1,0.9))
#'
#' # Now a data frame:
#' tst <- data.frame(mu = rnorm(1e4, 4, 1), sigma = rlnorm(1e4))
#' hdi(tst, 0.8)
#' apply(tst, 2, quantile, c(0.1,0.9))
#'
#' # For a function:
#' hdi(qgamma, 0.8, shape=2.5, rate=2)
#' # and the symmetrical 80% CrI:
#' qgamma(c(0.1, 0.9), 2.5, 2)
#'
#' # A severely bimodal distribution:
#' tst2 <- c(rnorm(1e5), rnorm(5e4, 7))
#' hist(tst2, freq=FALSE)
#' (hdiMC <- hdi(tst2))
#' segments(hdiMC[1], 0, hdiMC[2], 0, lwd=3, col='red')
#' # This is a valid 95% CrI, but not a Highest Density Interval
#'
#' dens2 <- density(tst2)
#' lines(dens2, lwd=2, col='blue')
#' (hdiD1 <- hdi(dens2))  # default allowSplit = FALSE; note the warning
#' segments(hdiD1[1], 0.01, hdiD1[2], 0.01, lty=3, col='blue')
#' # This is a valid 95% CrI, but not an HDI
#' (hdiD2 <- hdi(dens2, allowSplit=TRUE))
#' (ht <- attr(hdiD2, "height"))
#' segments(hdiD2[, 1], ht, hdiD2[, 2], ht, lwd=3, col='blue')
#' # This is the correct 95% HDI.
#' @export
hdi <- function(object, credMass=0.95, ...) UseMethod("hdi")

#' @describeIn hdi Calculates the hdi given a numeric vector
#' @export
hdi.default <- function(object, credMass=0.95, ...) {
  if(!is.numeric(object))
    stop(paste("No applicable method for class", class(object)))
  if(is.na(credMass) || length(credMass) != 1 || credMass <= 0 || credMass >= 1)
    stop("credMass must be in 0 < credMass < 1")
  if(all(is.na(object)))
    return(c(lower = NA_real_, upper = NA_real_))
  # This is Mike's code from way back:
  x <- sort(object)  # also removes NAs
  n <- length(x)
  # exclude <- ceiling(n * (1 - credMass)) # Not always the same as...
  exclude <- n - floor(n * credMass)       # Number of values to exclude
  low.poss <- x[1:exclude]             # Possible lower limits...
  upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
  best <- which.min(upp.poss - low.poss)      # Combination giving the narrowest interval
  result <- c(lower = low.poss[best], upper = upp.poss[best])

  attr(result, "credMass") <- credMass
  return(result)
}

#' @describeIn hdi Calculates the hdi of each column in a numeric matrix
#' @export
hdi.matrix <- function(object, credMass=0.95, ...) {
  result <- apply(object, 2, hdi.default, credMass=credMass, ...)
  attr(result, "credMass") <- credMass
  return(result)
}

#' @describeIn hdi Calculates the hdi of each column in a data frame
#' @export
hdi.data.frame <- function(object, credMass=0.95, ...)
  hdi.matrix(as.matrix(object), credMass=credMass, ...)

#' @describeIn hdi Calculates the hdi of each parameter in an mcmc.list object
#' @export
hdi.mcmc.list <- function(object, credMass=0.95, ...)
  hdi.matrix(as.matrix(object), credMass=credMass, ...)

#' @describeIn hdi Calculates the hdi given a inverse cumulative density function
#' @export
hdi.function <- function(object, credMass=0.95, tol, ...)  {
  if(is.na(credMass) || length(credMass) != 1 || credMass <= 0 || credMass >= 1)
    stop("credMass must be in 0 < credMass < 1")
  if(missing(tol))
    tol <- 1e-8
  if(class(try(object(0.5, ...), TRUE)) == "try-error")
    stop(paste("Incorrect arguments for the inverse cumulative density function",
        substitute(object)))
  # cf. code in Kruschke 2011 p630
   intervalWidth <- function( lowTailPr , ICDF , credMass , ... ) {
      ICDF( credMass + lowTailPr , ... ) - ICDF( lowTailPr , ... )
   }
   optInfo <- optimize( intervalWidth , c( 0 , 1.0 - credMass) , ICDF=object ,
                        credMass=credMass , tol=tol , ... )
   HDIlowTailPr <- optInfo$minimum
   result <- c(lower = object( HDIlowTailPr , ... ) ,
	            upper = object( credMass + HDIlowTailPr , ... ) )
  attr(result, "credMass") <- credMass
  return(result)
}

#' @describeIn hdi Calculates the hdi given the output from a call to \code{density()}
#' @export
hdi.density <- function(object, credMass=0.95, allowSplit=FALSE, ...) {
  if(is.na(credMass) || length(credMass) != 1 || credMass <= 0 || credMass >= 1)
    stop("credMass must be in 0 < credMass < 1")
  sorted = sort( object$y , decreasing=TRUE )
  heightIdx = min( which( cumsum( sorted) >= sum(object$y) * credMass ) )
  height = sorted[heightIdx]
  indices = which( object$y >= height )
  # HDImass = sum( object$y[indices] ) / sum(object$y)
  gaps <- which(diff(indices) > 1)
  if(length(gaps) > 0 && !allowSplit) {
    # In this case, return shortest 95% CrI
    warning("The HDI is discontinuous but allowSplit = FALSE;
    the result is a valid CrI but not HDI.")
    cumul <- cumsum(object$y) / sum(object$y)
    upp.poss <- low.poss <- which(cumul < 1 - credMass)
    for (i in low.poss)
      upp.poss[i] <- min(which(cumul > cumul[i] + credMass))
    # all(cumul[upp.poss] - cumul[low.poss] > credMass) # check
    width <- upp.poss - low.poss
    best <- which(width == min(width))  # usually > 1 value due to ties
    result <- c(lower = mean(object$x[low.poss[best]]),
                upper = mean(object$x[upp.poss[best]]))
  } else {
    begs <- indices[c(1, gaps + 1)]
    ends <- indices[c(gaps, length(indices))]
    result <- cbind(begin = object$x[begs], end = object$x[ends])
    if(!allowSplit)
      names(result) <- c("lower", "upper")
  }
  attr(result, "credMass") <- credMass
  attr(result, "height") <- height
  return(result)
}





