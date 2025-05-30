% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zetadiv.R
\name{Zeta.order.ex}
\alias{Zeta.order.ex}
\title{Expectation of zeta diversity for a specific number of assemblages or sites}
\usage{
Zeta.order.ex(
  data.spec,
  order = 1,
  sd.correct = TRUE,
  rescale = FALSE,
  empty.row = "empty"
)
}
\arguments{
\item{data.spec}{Site-by-species presence-absence data frame, with sites as rows and species as columns.}

\item{order}{Specific number of assemblages or sites at which zeta diversity is computed.}

\item{sd.correct}{Boolean value (TRUE or FALSE) indicating if the standard deviation must be computed with an unbiased estimator (using the number of site combinations - 1 as the denominator) or not (using the number of site combinations as the denominator).}

\item{rescale}{Boolean value (TRUE or FALSE) indicating if the zeta values should be divided by \eqn{\zeta_1}, to get a range of values between 0 and 1.}

\item{empty.row}{Determines how to handle empty rows, i.e. sites with no species. Such sites can cause underestimations of zeta diversity. Options are  "empty" to let the data untreated or "remove" to remove the empty rows.}
}
\value{
\code{zeta.order.ex}  returns a list containing the following components:

\item{zeta.order}{The number of assemblages or sites for which the zeta diversity was computed.}

\item{combinations}{The number of possible combinations of sites for the chosen order.}

\item{zeta.val}{The zeta diversity values.}

\item{zeta.val.sd}{The standard deviation of zeta diversity.}
}
\description{
Computes the expectation of zeta diversity, the number of species shared by multiple assemblages, for a specific order (number of assemblages or sites) using a formula based on the occupancy of the species.
}
\details{
\code{Zeta.order.ex} is much faster than \code{\link{Zeta.order.mc}} to compute the exact value of zeta diversity when the number of species is lower than \eqn{C^N_{i}}, where \emph{N} is the total number of sites and \emph{i} is the order of zeta.

\code{sd.correct} should be set to \code{TRUE} if the assemblages represent a subsample of the whole system. It can be set to \code{FALSE} if the sampling is exhaustive, for example in case of a continuous regular grid covering the whole study area.
}
\examples{

utils::data(bird.spec.coarse)
data.spec.bird <- bird.spec.coarse[,3:193]

zeta.bird <- Zeta.order.ex(data.spec.bird, order = 3)
zeta.bird

##########

utils::data(Marion.species)
data.spec.marion <- Marion.species[,3:33]

zeta.marion <- Zeta.order.ex(data.spec.marion, order = 3)
zeta.marion

}
\references{
Hui C. & McGeoch M.A. (2014). zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.

McGeoch M. A., Latombe G., Andrew N. R., Nakagawa S., Nipperess D. A., Roige M., Marzinelli E. M., Campbell A. H., Verges A., Thomas T., Steinberg P. D., Selwood K. E., Henriksen M. V. & Hui C. (2019). Measuring continuous compositional change using decline and decay in zeta diversity. \emph{Ecology}, 100(11), e02832.
}
\seealso{
\code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.decline.mc}}
}
