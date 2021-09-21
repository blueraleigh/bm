#' Phylogenetic independent contrasts
#'
#' Performs phylogenetic independent contrasts under a Brownian motion
#' process model.
#'
#' @param x A named vector of character state data having type \code{numeric}.
#' @param phy An object of class \code{tree}.
#' @return A list with two components:
#' \describe{
#' \item{pic}{A matrix of phylogenetic independent contrasts. Each row contains
#' the contrast value (first column) and contrast variance (second column) for
#' a single internal node. Row names correspond to internal node indices.}
#' \item{ace}{A vector of ancestral character state estimates for internal
#' nodes. Names correspond to internal node indices. Note that these are not
#' maximum likelihood estimates.}
#' }
#' Note that the returned list has a \code{rate} attribute that contains the
#' maximum likelihood estimate of the Brownian motion rate parameter.
#' @example
#' data(squamatatree)
#' data(squamatamass)
#' phy = read.newick(text=squamatatree)
#' fit = bm.pic(squamatamass, phy)
#' # log likelihood of data under Brownian motion
#' sum(dnorm(fit[[1]][,1], sd=sqrt(attr(fit[[1]], 'rate')*fit[[1]][,2]), log=TRUE))
#' # ancestral state estimates
#' fit[[2]]
bm.pic = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.numeric(x))

    if (length(setdiff(tiplabels(phy), names(x))))
        stop("Some terminal nodes are missing data")

    x = x[tiplabels(phy)]

    ans = .Call(C_bm_pic, x, phy)

    idx = as.character(root(phy):Nnode(phy))

    dimnames(ans[[1]]) = list(idx, c("u", "u.var"))
    names(ans[[2]]) = idx

    names(ans) = c("pic", "ace")

    return (ans)
}
