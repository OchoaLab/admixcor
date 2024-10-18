# identify the individuals that are best candidates to be in vertices
# most of the code is for the edge cases of ties, which are disambiguated with random noise.  It's very important that indexes be unique to avoid direct colinearity in a matrix inversion later
vertex_inds <- function ( Q ) {
    # number of ancestries
    k_subpops <- ncol( Q )
    
    # this simple code finds the individual with the most of each ancestry
    indexes <- apply( Q, 2, which.max )
    
    # see if there's repeated individuals, and fix that case (in a loop involving random steps)
    while ( length( unique( indexes ) ) != k_subpops ) {
        # actually find the ones that are repeated, to add noise to those rows only
        x <- table( indexes )
        indexes_ties <- x[ x > 1 ]
        # handle the case that there's more than one of these, extremely unlikely but meh
        n_ties <- length( indexes_ties )
        # to break ties, just add very small amounts of noise at those positions
        noise <- matrix(
            stats::rnorm( k_subpops * n_ties, sd = 1e-6 ),
            nrow = n_ties,
            ncol = k_subpops
        )
        # hack works if noise has sample mean zero
        noise <- noise - rowMeans( noise )
        # create a copy to add the noise to
        Q2 <- Q
        Q2[ indexes_ties, ] <- Q2[ indexes_ties, ] + noise 
        # recalculate, repeat until there's no ties
        indexes <- apply( Q2, 2, which.max )
    }
    return( indexes )
}
