# goal is to align an estimated matrix Q2 to a reference or true matrix Q
# this approach is based on:
# - similarity metric: RMSD
# - algorithm: brute force enumeration of all permutations! (uses gtools package)
#   - ok performance at k=10 when fast=TRUE
#   - fast=FALSE is most direct/naive version, for internal testing only (scales extremely poorly for K=10)
# originally by Amika Sood, modified by Alex Ochoa
#
# tests with R01 data, approach ADM, all viiiaR5, with some lines uncommented
# Rscript sim-bed-01-align-Q.R -a ADM-$k -k $k
#  k    orig   fast  fast2  heur (*=RMSD gap 0)
#  5   0.011  0.001  0.002 0.028*
#  6   0.046  0.007  0.008 0.028*
#  7   0.271  0.050  0.061 0.029*
#  8   2.397  0.391  0.490 0.028 (gap 0.0017038)
#  9  22.929  3.579  4.511 0.028 (gap 0.0007097)
# 10 229.040 36.267 46.199 0.029 (gap 0.0099809)

#' @export
align_Q <- function( Q, Q2, verbose = FALSE, fast = TRUE, fast2 = FALSE ) {
    # number of ancestral populations
    k <- ncol(Q)
    # number of individuals
    n <- nrow(Q)
    
    # validations
    if ( ncol(Q2) != k )
        stop( 'Q and Q2 dimensions disagree! K = ', k, ' vs ', ncol(Q2) )
    if ( nrow(Q2) != n )
        stop( 'Q and Q2 dimensions disagree! n = ', n, ' vs ', nrow(Q2) )
    
    # calculate RMSD contributions pairwise matrix
    # ought to be relatively fast
    if ( fast )
        rmsd_mat <- rmsd_Q_mat( Q, Q2 )
    
    # enumerate all permutations
    if (verbose) message( 'permutations' )
    permut <- gtools::permutations(n=k,r=k)
    err <- NA
    i_max <- factorial(k)
    for(i in 1:i_max){
        # current permutation
        v <- permut[ i, ]
        if ( fast ) {
            # these versions factor out the computation of distances between columns, so they're fastest by a lot
            if ( fast2 ) {
                # avoiding the loop, but is slower by having to permute entire matrices
                e <- sum( diag( rmsd_mat[ , v ] ) )
            } else {
                # fastest version: direct sum of desired elements only
                # (have to walk [ j, v[j] ] path in explicit loop)
                e <- 0
                for ( j in 1 : k ) {
                    e <- e + rmsd_mat[ j, v[ j ] ]
                }
            }
        } else {
            # direct RMSD for permuted matrix (super slow!)
            e <- rmsd_Q( Q, Q2[ , v ] )
        }
        # remember best
        if ( is.na( err ) || err > e ) {
            err <- e
            # Q2b <- Q2p # not returned!
            indexes <- v
        }
    }
    
    # return permutation indexes only, to reorder other things like P and Psi outside
    return( indexes )
}
