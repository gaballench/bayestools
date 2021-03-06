% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/measureSensit.R
\name{measureSensit}
\alias{measureSensit}
\title{Calculate the Intersection Between Two Densities}
\usage{
measureSensit(d1, d2, splits = 500, rawData = c(TRUE, TRUE), plot = TRUE,
  x_limit = "auto", colors = c("red", "blue", "gray"), ...)
}
\arguments{
\item{d1, d2}{Either two vectors of empirical (i.e., MCMC-produced) values OR a \code{data.frame}/\code{matrix} with columns x and y for values fitted to a density from which to calculate areas. If \code{rawData} is set to \code{TRUE} in any instance, the data must be placed in vectors and not multidimensional objects.}

\item{splits}{A numerical argument controling the number of subdivisions of the intersection area for numerical integration}

\item{rawData}{Are d1 and/or d2 raw data for which a density should be calculated? A vector of length two containing logical values indicating whenther any of the arguments d1 or d2 are raw data or whether the user is inputing already calculated densities (e.g., the output from the density, curve, or dDIST functions, or any two-dimension object with x and y values)}

\item{plot}{Should a plot be produced?}

\item{colors}{A vector of three colors, namely, color of the \code{d1} density (e.g., the prior), color of the \code{d2} density e.g., the posterior), and color of the intersection.}

\item{...}{Further arguments to pass to the graphical functions such as \code{lines} and \code{plot} internally (e.g., \code{main}, \code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, etc.).}
}
\value{
A numeric vector with the value of the intersection between two densities. As a side effect, a plot is produced to an active (or new) graphical device.
}
\description{
Calculate the Intersection Between Two Densities
}
\details{
Sensitivity is measured as the overlapping portion between two densities following Ballen (in prep). It has a value between 0 and 1. The values of the vector rawData determine the behavior of the function and therefore attention must be paid to their consistence with the nature of arguments d1 and d2. Despite the function was designed in order to allow to quantify posterior sensitivity to the prior, this can be used to quantify any overlap between two given densities and for any other purpose.
}
\examples{
\dontrun{
# Set seed and colors to use in plots in the order: Prior, posterior, and intersection
set.seed(1985)
colors <- c("red", "blue", "lightgray")
# Sensitivity in two identical distributions
below <- measureSensit(d1 = rnorm(1000000, mean = 0, 1),
                       d2 = rnorm(1000000, mean = 0, 1),
                       main = "Comp. dependence",
                       colors = colors)
legend(x = "topright", legend = round(below, digits = 2))
# Sensitivity in two distributions partially overlapping
below <- measureSensit(d1 = rnorm(1000000, mean = 3, 1),
                       d2 = rnorm(1000000, mean = 0, 1),
                       main = "Partial dependence",
                       colors = colors)
legend(x = "topright", legend = round(below, digits = 2))
# Sensitivity in two completely-different distributions
below <- measureSensit(d1 = rnorm(1000000, mean = 8, 1),
                       d2 = rnorm(1000000, mean = 0, 1),
                       main = "Comp. independence",
                       colors = colors)
legend(x = "topright", legend = round(below, digits = 2))
# Don't plot, just return the intersection
measureSensit(d1 = rnorm(1000000, mean = 3, 1),
              d2 = rnorm(1000000, mean = 0, 1),
              plot = FALSE)
}
}
