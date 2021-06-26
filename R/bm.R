#' @example
#' data(squamatatree)
#' data(squamatamass)
#' phy = read.newick(text=squamatatree)
#' phy = drop.tip(phy, setdiff(names(squamatamass), tiplabels(phy)))
#' x = bm.pic(squamatamass, phy)
#' # log likelihood of data under Brownian motion
#' sum(dnorm(x[[1]][,1], sd=sqrt(attr(x[[1]], 'rate')*x[[1]][,2]), log=TRUE))
#' # ancestral state estimates
#' x[[2]]
bm.pic = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.numeric(x))
    x = x[tiplabels(phy)]

    if (length(setdiff(tiplabels(phy), names(x))))
        stop("Some terminal nodes are missing data")

    ans = .Call(C_bm_pic, x, phy)

    idx = as.character(root(phy):Nnode(phy))

    dimnames(ans[[1]]) = list(idx, c("u", "u.var"))
    names(ans[[2]]) = idx

    names(ans) = c("pic", "ace")

    return (ans)
}
