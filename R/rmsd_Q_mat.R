# compute an RMSD-related distance matrix between all columns of two Q matrices, to help find best permutation
# each value is not an RMSD measure, but a per-pair contribution which is an RMSD only after normalizing, adding across all pairs, and taking the square root
# INTERNAL, only for use in align_Q* versions above
rmsd_Q_mat <- function( Q, Q2 ) {
    # number of ancestral populations
    k <- ncol(Q)
    # initialize RMSD matrix
    rmsd_mat <- matrix( 0, ncol = k, nrow = k )
    # navigate all pairs
    # note that matrix is asymetric
    for ( i in 1 : k ) {
        # extract column
        qi <- Q[ , i ]
        for ( j in 1 : k ) {
            # extract column
            q2j <- Q2[ , j ]
            # calculate and save value for Q2's column
            # NOTE: didn't use "norm()" because it requires matrix inputs
            rmsd_mat[ i, j ] <- sum( ( qi - q2j )^2 )
        }
    }
    # to make it more RMSD like, could normalize by n and k, but meh, minimum doesn't change either way
    # return this matrix of values
    return( rmsd_mat )
}
