% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bm.R
\name{bm.pic}
\alias{bm.pic}
\title{Phylogenetic independent contrasts}
\usage{
bm.pic(x, phy)
}
\arguments{
\item{x}{A named vector of character state data having type \code{numeric}.}

\item{phy}{An object of class \code{tree}.}
}
\value{
A list with two components:
\describe{
\item{pic}{A matrix of phylogenetic independent contrasts. Each row contains
the contrast value (first column) and contrast variance (second column) for
a single internal node. Row names correspond to internal node indices.}
\item{ace}{A vector of ancestral character state estimates for internal
nodes. Names correspond to internal node indices. Note that these are not
maximum likelihood estimates.}
}
Note that the returned list has a \code{rate} attribute that contains the
maximum likelihood estimate of the Brownian motion rate parameter.
}
\description{
Performs phylogenetic independent contrasts under a Brownian motion
process model.
}
\examples{
data(squamatatree)
data(squamatamass)
phy = read.newick(text=squamatatree)
fit = bm.pic(squamatamass, phy)
# log likelihood of data under Brownian motion
sum(dnorm(fit[[1]][,1], sd=sqrt(attr(fit[[1]], 'rate')*fit[[1]][,2]), log=TRUE))
# ancestral state estimates
fit[[2]]
}
