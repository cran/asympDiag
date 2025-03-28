% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/envelope.R
\name{plot.AD_envelope}
\alias{plot.AD_envelope}
\title{Envelope Plot}
\usage{
\method{plot}{AD_envelope}(
  x,
  colors = getOption("asympDiag.plot.AD_envelope.colors"),
  xlab = "Expected quantiles",
  ylab = "Observed quantiles",
  distribution = function(p) stats::qnorm((1 + p)/2),
  ylim = base::range(c(x$observed, x$lower, x$upper), na.rm = TRUE, finite = TRUE),
  ...
)
}
\arguments{
\item{x}{AD_envelope object, usually the result of \code{\link[=envelope]{envelope()}}}

\item{colors}{Vector of length 2, with color for points outside and inside
the envelope band, respectively.}

\item{xlab}{The label for the x-axis.}

\item{ylab}{The label for the y-axis.}

\item{distribution}{quantile function for reference theoretical distribution.}

\item{ylim}{the y limits of the plot.}

\item{...}{extra arguments passed to \link[graphics:plot.default]{graphics::plot}}
}
\value{
No return value, called for side effects
}
\description{
Plot AD_envelope
}
\details{
Create envelope plot. The expected quantile, by default, is from a half-normal
distribution, as the default residuals from \code{\link[=envelope_residual]{envelope_residual()}} are absolute.
You may replace it with the \code{distribution} argument.
}
