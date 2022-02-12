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
#' @examples
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


#' Compute tip rates using phylogenetic independent contrasts
#'
#' As the number of terminal nodes in the phylogeny grows large the
#' average tip rate converges on the maximum likelihood estimate of the
#' Brownian motion variance parameter.
#'
#' @param x A named vector of character state data having type \code{numeric}.
#' @param phy An object of class \code{tree}.
#' @return A vector of tip rates
#' @examples
#' data(squamatatree)
#' data(squamatamass)
#' phy = read.newick(text=squamatatree)
#' x = bm.tiprate(squamatamass, phy)
#' # average tiprate
#' mean(x) # 0.02721007
#' # maximum likelihood rate estimate
#' attr(bm.pic(squamatamass, phy)[[1]], "rate") # 0.02721699
bm.tiprate = function(x, phy) {
    n = Ntip(phy)
    y = numeric(n)
    pic = bm.pic(x, phy)
    for (i in 1:n) {
        anc = ancestors(phy)[[i]]
        for (j in 1:length(anc)) {
            cont = pic[[1]][as.character(anc[j]), ]
            # squared contrast scaled by contrast variance
            d = (cont[1]*cont[1]) / cont[2]
            y[i] = y[i] + d / 2^j
        }
    }
    y
}


#' Compute multivariate tip rates using phylogenetic independent contrasts
#'
#' As the number of terminal nodes in the phylogeny grows large the
#' average tip rate converges on the maximum likelihood estimate of the
#' Brownian motion covariance matrix
#'
#' @param x A matrix (with row names) of character state data having type \code{numeric}.
#' @param phy An object of class \code{tree}.
#' @return A list of matrix-valued tip rates
bm.mvtiprate = function(x, phy) {
    p = ncol(x)
    n = Ntip(phy)
    U = apply(x, 2, bm.pic, phy=phy)
    U = do.call(cbind, lapply(U, function(u) u[[1]][,1] / sqrt(u[[1]][,2])))
    y = vector("list", n)
    for (i in 1:n) {
        y[[i]] = matrix(0, p, p)
        anc = ancestors(phy)[[i]]
        for (j in 1:length(anc)) {
            u = U[as.character(anc[j]), ]
            d = u %*% t(u)
            y[[i]] = y[[i]] + d / 2^j
        }
    }
    y
}
