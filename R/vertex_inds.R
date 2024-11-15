# identify the individuals that are best candidates to be in vertices
# most of the code is for the edge cases of ties, there are two startegies for select some vertices to keep (or not) and the rest (or all) are returned as NAs (no vertex assigned; outside these ancestries are identified to not stretch/modify at all).  It's very important that indexes be unique to avoid direct colinearity in a matrix inversion later
vertex_inds <- function ( Q, ties_none = FALSE ) {
    # number of ancestries
    k_subpops <- ncol( Q )
    
    # this simple code finds the individual with the most of each ancestry
    indexes <- apply( Q, 2, which.max )

    # see if there's repeated individuals, and fix that case (in a loop involving random steps)
    if ( length( unique( indexes ) ) != k_subpops ) {
        # actually find the ones that are repeated
        x <- table( indexes )
        indexes_ties <- as.integer( names( x[ x > 1 ] ) )
        
        if ( ties_none ) {
            # one solution is to keep none of these, return all as NAs
            indexes[ indexes %in% indexes_ties ] <- NA
        } else {
            # alternatively, decide which of the two or more ancestries to keep this individual assigned to, based on majority, set the rest to NA
            # navigate individual that appeared for multiple ancestries
            for ( index in indexes_ties ) {
                # these were the ancestries that claimed this individual
                ks <- which( indexes == index )
                # this is the most common ancestry, between the competitors, for this individual
                # (turns out if we look at all ancestries, another one (not in `ks`) can have even higher values!)
                k_best <- ks[ which.max( Q[ index, ks ] ) ]
                # in this case remove k_best from ks, the rest get turned into NAs
                ks <- setdiff( ks, k_best )
                indexes[ ks ] <- NA
            }
        }
    }
    return( indexes )
}
