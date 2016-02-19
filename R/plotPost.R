#' Graphic display of a posterior probability distribution
#'
#' Plot the posterior probability distribution for a single parameter from a
#' vector of samples, typically from an MCMC process, with appropriate summary
#' statistics.
#'
#' The data are plotted either as a histogram (above) or, if \code{showCurve =
#' TRUE}, as a fitted kernel density curve (below). Either the mean or the mode
#' of the distribution is displayed, depending on the parameter \code{showMode.}
#' The Highest Density Interval (HDI) is shown as a horizontal bar, with labels
#' for the ends of the interval.
#'
#' \if{html}{\figure{plotPost1.jpg} }
#' \if{latex}{\figure{plotPost1.jpg}{options: width=5cm}}
#' \cr
#' \cr
#' \if{html}{\figure{plotPost2.jpg} }
#' \if{latex}{\figure{plotPost2.jpg}{options: width=5cm}}
#'
#' If values for a ROPE are supplied, these are shown as dark red vertical
#' dashed lines, together with the percentage of probability mass within the
#' ROPE. If a comparison value (\code{compVal}) is supplied, this is shown as a
#' vertical green dotted line, together with the probability mass below and
#' above this value.
#'
#' @param paramSampleVec A vector of samples drawn from the target distribution.
#' @param credMass the probability mass to include in credible intervals, or
#'   NULL to suppress plotting of credible intervals.
#' @param compVal a value for comparison with those plotted.
#' @param ROPE a two element vector, such as \code{c(-1, 1)}, specifying the
#'   limits of the Region Of Practical Equivalence.
#' @param HDItextPlace a value in [0,1] that controls the horizontal position of
#'   the labels at the ends of the HDI bar.
#' @param showMode logical: if TRUE, the mode is displayed instead of the mean.
#' @param showCurve logical: if TRUE, the posterior density will be represented
#'   by a kernel density function instead of a histogram.
#' @param \dots graphical parameters and the \code{breaks} parameter for the
#'   histogram.
#' @return Returns an object of class \code{histogram} invisibly. Used for its
#'   plotting side-effect.
#'
#' @note The origin of this function is
#'   \href{http://cran.r-project.org/package=BEST}{the BEST
#'   package} which is based on Kruschke(2015, 2013).
#'
#' @author John Kruschke, modified by Mike Meredith
#' @seealso For details of the HDI calculation, see \code{\link{hdi}}.
#' @examples
#' # Generate some data
#' tst <- rnorm(1e5, 3, 1)
#' plotPost(tst)
#' plotPost(tst, col='wheat', border='magenta')
#' plotPost(tst, credMass=0.8, ROPE=c(-1,1), xlab="Response variable")
#' plotPost(tst, showMode=TRUE, showCurve=TRUE, compVal=5.5)
#'
#' # For integers:
#' tst <- rpois(1e5, 12)
#' plotPost(tst)
#'
#' # A severely bimodal distribution:
#' tst2 <- c(rnorm(1e5), rnorm(5e4, 7))
#' plotPost(tst2)                  # A valid 95% CrI, but not HDI
#' plotPost(tst2, showCurve=TRUE)  # Correct 95% HDI
#' @references
#' Kruschke, J. K. (2015) \emph{Doing Bayesian data analysis, second
#'   edition: A tutorial with R, JAGS, and Stan.} Waltham, MA: Academic Press /
#'   Elsevier.
#'
#' Kruschke, J. K. (2013) Bayesian estimation supersedes the t test.
#'   \emph{Journal of Experimental Psychology: General}, \bold{142(2)}, 573.
#' @export
plotPost <-
  function(paramSampleVec, credMass = 0.95, compVal = NULL, ROPE = NULL,
           HDItextPlace = 0.7, showMode = FALSE, showCurve = FALSE, ...) {

  # Does a plot for a single parameter. Called by plot.BEST but also exported.
  # Returns a histogram object invisibly.
  # This stuff should be in the ... argument:
  #   yaxt="n", ylab="", xlab="Parameter", main="", cex.lab=1.5, cex=1.4,
  #   xlim=range(compVal, paramSampleVec), col="skyblue", border="white",
  #   breaks=NULL

  # Deal with ... argument:
  dots <- list(...)
  if(length(dots) == 1 && class(dots[[1]]) == "list")
    dots <- dots[[1]]
  defaultArgs <- list(xlab=deparse(substitute(paramSampleVec)),
    yaxt="n", ylab="", main="", cex.lab=1.5,
    cex=1.4, col="skyblue", border="white", bty="n", lwd=5, freq=FALSE,
    xlim=range(compVal, HDInterval::hdi(paramSampleVec, 0.99)))
  useArgs <- modifyList(defaultArgs, dots)

  # Get breaks argument
  breaks <- dots$breaks
  if (is.null(breaks)) {
    if (all(paramSampleVec == round(paramSampleVec))) { # all integers
      breaks <- seq(min(paramSampleVec), max(paramSampleVec) + 1) - 0.5
    } else {
      by <- diff(HDInterval::hdi(paramSampleVec))/18
      breaks <- unique(c( seq( from=min(paramSampleVec), to=max(paramSampleVec),
                     by=by), max(paramSampleVec) ))
    }
  }
  histinfo <- hist(paramSampleVec, breaks=breaks, plot=FALSE)
  histinfo$xname <- useArgs$xlab

  oldpar <- par(xpd=TRUE) ; on.exit(par(oldpar))

  if (showCurve) {
    densCurve <- density( paramSampleVec, adjust=2 )
    selPlot <- names(useArgs) %in%
      c(names(as.list(args(plot.default))), names(par(no.readonly=TRUE)))
    plotArgs <- useArgs[selPlot]
    plotArgs$x <- densCurve$x
    plotArgs$y <- densCurve$y
    plotArgs$type <- "l"
    plotArgs$xpd <- FALSE
    do.call(plot, plotArgs, quote=TRUE)
    abline(h=0, col='grey', xpd=FALSE)
    # Display the HDI.
    if(!is.null(credMass)) {
      HDI <- HDInterval::hdi(densCurve, credMass, allowSplit=TRUE)
      ht <- attr(HDI, "height")
      segments(HDI[, 1], ht, HDI[, 2], ht, lwd=4, lend='butt')
      segments(HDI, 0, HDI, ht, lty=2)
      text( mean(HDI), ht, bquote(.(100*credMass) * "% HDI" ),
            adj=c(.5,-1.7), cex=useArgs$cex )
      text( HDI, ht, bquote(.(signif(HDI, 3))),
            pos=3, cex=useArgs$cex )
    }
  } else {
    plot.histogram.args.names <- c("freq", "density", "angle", "border",
      "main", "sub", "xlab", "ylab", "xlim", "ylim", "axes", "labels",
      "add") # plot.histogram not exported, so need to cheat!
    selPlot <- names(useArgs) %in%
      c(plot.histogram.args.names, names(par(no.readonly=TRUE)))
    plotArgs <- useArgs[selPlot]
    plotArgs$lwd <- 1
    plotArgs$x <- histinfo
    do.call(plot, plotArgs, quote=TRUE)
    # Display the HDI.
    if(!is.null(credMass)) {
      HDI <- HDInterval::hdi( paramSampleVec, credMass )
      lines(HDI, c(0,0), lwd=4, lend='butt')
      text( mean(HDI), 0, bquote(.(100*credMass) * "% HDI" ),
            adj=c(.5,-1.7), cex=useArgs$cex )
      text( HDI[1], 0, bquote(.(signif(HDI[1],3))),
            adj=c(HDItextPlace,-0.5), cex=useArgs$cex )
      text( HDI[2], 0, bquote(.(signif(HDI[2],3))),
            adj=c(1.0-HDItextPlace,-0.5), cex=useArgs$cex )
    }
  }


  # Display mean or mode:
  cenTendHt <- 0.9 * max(histinfo$density)
  if ( showMode==FALSE ) {
      meanParam <- mean( paramSampleVec )
      text( meanParam, cenTendHt,
            bquote(mean==.(signif(meanParam,3))), adj=c(.5,0), cex=useArgs$cex )
  } else {
      dres <- density( paramSampleVec )
      modeParam <- dres$x[which.max(dres$y)]
      text( modeParam, cenTendHt,
            bquote(mode==.(signif(modeParam,3))), adj=c(.5,0), cex=useArgs$cex )
  }
  # Display the comparison value.
  if ( !is.null( compVal ) ) {
    cvHt <- 0.7 * max(histinfo$density)
    cvCol <- "darkgreen"
    pcgtCompVal <- round( 100 * sum( paramSampleVec > compVal )
                          / length( paramSampleVec ) , 1 )
     pcltCompVal <- 100 - pcgtCompVal
     lines( c(compVal,compVal), c(0.96*cvHt,0),
            lty="dotted", lwd=1, col=cvCol )
     text( compVal, cvHt,
           bquote( .(pcltCompVal)*"% < " *
                   .(signif(compVal,3)) * " < "*.(pcgtCompVal)*"%" ),
           adj=c(pcltCompVal/100,0), cex=0.8*useArgs$cex, col=cvCol )
  }
  # Display the ROPE.
  if ( !is.null( ROPE ) ) {
    ROPEtextHt <- 0.55 * max(histinfo$density)
    ropeCol <- "darkred"
     pcInROPE <- ( sum( paramSampleVec > ROPE[1] & paramSampleVec < ROPE[2] )
                          / length( paramSampleVec ) )
     lines( c(ROPE[1],ROPE[1]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=2,
            col=ropeCol )
     lines( c(ROPE[2],ROPE[2]), c(0.96*ROPEtextHt,0), lty="dotted", lwd=2,
            col=ropeCol)
     text( mean(ROPE), ROPEtextHt,
           bquote( .(round(100*pcInROPE))*"% in ROPE" ),
           adj=c(.5,0), cex=1, col=ropeCol )
  }

  return(invisible(histinfo))
}
