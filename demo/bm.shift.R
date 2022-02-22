### Fitting multi-rate Brownian motion models of character evolution
### to empirical data using bm::bm.shift

library(bm)

# Data on turtle body size from
# Eastman, J.M., Alfara, M.E., Joyce, P., Hipp, A.L., and L.J. Harmon. 2011. 
# Evolution 65(12): 3578-3589
data(chelonia)

x = chelonia$dat
phy = phylo::read.newick(text=ape::write.tree(chelonia$phy))

# A named vector of class ‘numeric’ with log-transformed body size data for 
# 226 turtles
head(x)

# As a first exploration of the data we'll simply plot the distribution
# of contrasts
pic = bm::bm.pic(x, phy)

# The first component of this list is a matrix with independent contrasts
# stored in the first column and the contrast variances stored in the second/
# The rownames of this matrix correspond to the node indices to which the
# contrasts apply. For 'tree' class phylogenies used in bm mk package, node indices
# are identical to the node indices of an ape 'phylo' object in 'cladewise'
# node order
#
# The second component of the list are ancestral character state estimates.
# These are not maximum likelihood estimates, but they are locally parsimonious
# and made during the course of the contrast algorithm.


# We'll use these results to color branches according to their
# inferred mass, with warmer colors representing heavier masses.
node.state = c(x, pic$ace)
rate.bin = findInterval(node.state, seq(min(node.state), max(node.state), length.out=33))
pal = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", 
        "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
edge.color = colorRampPalette(pal)(33)[rate.bin]

# Note that the plot function invisibly returns the plotting coordinates,
# which we'll store for later annotations
L = plot(phy, edge.color=edge.color, lwd=1.5)
title("Contrast-inferred turtle body size evolution")

u = numeric(phylo::Nnode(phy))
u[phylo::root(phy):phylo::Nnode(phy)] = abs(pic$pic[, 1]) / sqrt(pic$pic[, 2])

# Now we want to highlight the nodes that show large (standardized) contrasts
points(L[[1]][, 1], L[[1]][, 3], cex=5*u, lwd=0.5)

# We can observe that large contrasts are definitely clustered in certain
# regions of the phylogeny that is suggestive of rate heterogeneity, so let's 
# fit a multi-rate model of body size evolution and see how it aligns with the
# distribution of contrasts
fit = bm::bm.shift(x, phy)

# This returns a list with three components.
#
# The first component is a vector of lineage-specific rates of character
# evolution that have been average over multiple possible rate-shift
# configurations. The second component is a set of weights assigned to
# different rate-shift configurations. These weights were used to compute
# the average rates in the first component. The final component is a function
# that can be used to retrieve the lineage-specific rates of character
# evolution associated with a particular rate-shift configuration.
str(fit)


# Let's first look at the model-averaged rates. Note that each value in
# this vector is the rate assigned to the edge leading to the node whose index 
# is the same as the vector index.
edge.rate = fit$avg.rates

# Let's discretize the rates and then map them to a set of colors, where
# warmer colors mean faster rates
rate.bin = findInterval(edge.rate, seq(0, max(edge.rate), length.out=33))
edge.color = colorRampPalette(
     c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))(33)[rate.bin]

plot(phy, edge.color=edge.color, lwd=1.5)
title("Model-averaged rates")
# Now we want to highlight the nodes that show large (standardized) contrasts
points(L[[1]][, 1], L[[1]][, 3], cex=5*u, col=1, lwd=0.5)

# And what you should see is that all the regions of phylogeny showing large 
# contrasts are inferred to have elevated rates of evolution compared to regions
# with small contrasts.


# We can repeat the exact same process using individual rate-shift configurations
# as well, not just the model-averaged rates. Let's look at the individual
# configuration that has the highest score (highest AIC weight)
best_index = which.max(fit$aic.weights)

# Note that we have to be a little careful here due to different indexing
# conventions between R and C. Rather than pass best_index directly, we need
# to pass best_index - 1. This is why the names associated with the
# fit$aic.weights vector are always one less than the vector position.
edge.rate = fit$rate(best_index-1L)

rate.bin = findInterval(edge.rate, seq(0, max(edge.rate), length.out=33))
edge.color = colorRampPalette(
     c("#BEBEBE","#00008B","#FF0000","#FFA500","#FFD700"))(33)[rate.bin]

# The overall pattern here is quite similar to the model-averaged pattern
plot(phy, edge.color=edge.color, lwd=1.5)
points(L[[1]][, 1], L[[1]][, 3], cex=5*u, col=1, lwd=0.5)
title("Best rate-shift configuration rates")
