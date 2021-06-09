#' @example
#' data(squamatatree)
#' data(squamatamass)
#' phy = read.newick(text=squamatatree)
#' phy = drop.tip(phy, setdiff(names(squamatamass), tiplabels(phy)))
#' x = mptp.bm(squamatamass, phy)
#' r8t = findInterval(log(x$avg.rates), seq(min(log(x$avg.rates)), max(log(x$avg.rates)), length.out=33))
#' edge.color = colorRampPalette(
#'    rev(c("#67001F",
#' "#B2182B",
#' "#D6604D",
#' "#F4A582",
#' "#FDDBC7",
#' "#F7F7F7",
#' "#D1E5F0",
#' "#92C5DE",
#' "#4393C3",
#' "#2166AC",
#' "#053061")))(33)[r8t]
#' plot(phy, edge.color=edge.color)
mptp.bm = function(x, phy) {
    stopifnot(!is.null(names(x)))
    stopifnot(is.numeric(x))
    x = x[tiplabels(phy)]
    obj = .Call(C_bm_shift, x, phy)

    list(
        avg.rates = obj[[1L]]
        , aic.weights = structure(obj[[2L]], names=0:(Nnode(phy)-1))
        , rate = function(index) {
            .Call(C_bm_shift_backtrack, as.integer(index), phy)
        }
    )
}
