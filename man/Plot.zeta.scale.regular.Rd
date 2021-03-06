% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zetadiv.R
\name{Plot.zeta.scale.regular}
\alias{Plot.zeta.scale.regular}
\title{Plotting of zeta diversity scaling with sample grain using hierarchical increases in grain size}
\usage{
Plot.zeta.scale.regular(
  zeta.scale.reg,
  size.init = 1,
  add = FALSE,
  ylim = NULL,
  col = "black"
)
}
\arguments{
\item{zeta.scale.reg}{A list generated by the function \code{Zeta.scale.regular}.}

\item{size.init}{initial Size of the plots before aggregation.}

\item{add}{Boolean value indicating if the graph must be plotted in a new graphics device or added to the active one.}

\item{ylim}{Numeric vectors of length 2, giving the range of y values.}

\item{col}{String indicating the color of the graph.}
}
\value{
A plot of the zeta diversity scaling with the mapping grain \code{n} x \code{n} (the number of sites combined to generate data at a coarser grain) on the x-axis and the value of zeta on the y-axis.
}
\description{
Plots the output of the function \code{Zeta.scale.regular}.
}
\examples{
utils::data(bird.spec.fine)
xy.bird <- bird.spec.fine[1:400,1:2]
data.spec.bird <- bird.spec.fine[1:400,3:192]

##sam = 25 is used here for fast execution, but a higher value is advised
zeta.scale.reg <- Zeta.scale.regular(xy.bird, data.spec.bird, n = 1:3, order = 3,
    sam = 25, normalize = "Jaccard",plot=FALSE)
dev.new()
Plot.zeta.scale.regular(zeta.scale.reg)
}
\references{
Hui C. & McGeoch M.A. (2014). Zeta diversity as a concept and metric that unifies incidence-based biodiversity patterns. \emph{The American Naturalist}, 184, 684-694.

Scheiner S.M., Chiarucci A., Fox G.A., Helmus M.R., McGlinn D.J. & Willig M.R. (2011). The underpinnings of the relationship of species richness with space and time. \emph{Ecological Monographs}, 81, 195-213.
}
\seealso{
\code{\link{Zeta.decline.mc}}, \code{\link{Zeta.order.mc}}, \code{\link{Zeta.decline.ex}}, \code{\link{Zeta.order.ex}},

\code{\link{Zeta.scale.regular}}, \code{\link{Zeta.scale.min.dist}}, \code{\link{rescale.regular}},

\code{\link{Plot.zeta.scale.min.dist}}
}
