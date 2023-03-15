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

# aligns to minimize RMSD as above, but this one uses a heuristic to be much faster for large K
# despite some limited 2nd-order search tricks, still doesn't always identify the best alignment in simple cases (i.e. k=5), so I'd more strongly favor the exact (exhaustive search) version `align_Q`
align_Q_heuristic <- function( Q, Q2 ) {
    # calculate RMSD contributions pairwise matrix
    # ought to be relatively fast
    rmsd_mat <- rmsd_Q_mat( Q, Q2 )
    
    # so the permutation ought to pick every column of Q2 exactly once (paired with each Q column, or row in rmsd_mat) and minimize the sum of the path through the RMSD matrix
    # so ideally there's a clear 1-1 alignment, and in case of ties or etc we can use the ranks to figure it out heuristically?
    # anyway, let's get the ordered preferences per column of Q
    k <- ncol(Q)
    # list is better, as we edit indexes each map will have different lengths
    rmsd_order <- vector( 'list', k )
    for ( i in 1 : k ) {
        rmsd_order[[ i ]] <- order( rmsd_mat[ i, ] )
    }

    # if the problem is easy, then the columns of Q and Q2 are already a 1-1 map (the first column of rmsd_order) and we're done!
    indexes <- sapply( rmsd_order, function(x) x[1] ) # get first element
    if ( length( indexes ) == length( unique( indexes ) ) )
        return( indexes )

    # else the obvious best match solution repeats some of the columns of Q2
    # the first part of the heuristic is to keep singletons fixed (the column in Q2 that was best match to a given column in Q and not a best match for any other columns of Q)
    counts <- table(indexes)
    indexes_done <- as.integer( names( counts[ counts == 1 ] ) )

    # clean up the done cases from order map
    # to think clearly as eliminating possibilities
    for ( i in 1 : k ) {
        rmsd_order_i <- rmsd_order[[ i ]]
        if ( rmsd_order_i[1] %in% indexes_done ) {
            # in this case the top choice is accepted, delete all other choices for visual clarity
            rmsd_order[[ i ]] <- rmsd_order_i[1]
        } else {
            # the top choice hasn't been accepted
            # however, remove fixed cases (from other columns) as options here
            rmsd_order[[ i ]] <- setdiff( rmsd_order_i, indexes_done )
        }
    }

    # now make decisions for the remaining clashes
    # this will be iterative
    while ( length(indexes_done) < k ) {
        # the heuristic fixes the clashing best-pair with the smallest score, forcing the others to unclash

        # first gather all the scores for clashing best pairs
        scores <- vector( 'integer', k )
        for ( i in 1 : k ) {
            # get the best pair index
            j <- rmsd_order[[ i ]][1]
            # do not score if j is already fixed
            if ( j %in% indexes_done ) {
                scores[i] <- NA
            } else {
                # this is the score we want to copy for this evaluation
                scores[i] <- rmsd_mat[i,j]
            }
        }

        # now we identify the minimum
        # this automatically excludes NAs, perfect because those are fixed already
        # in the unlikely case of ties, always returns first match (always scalar)
        i_fix <- which.min( scores )

        ### TODO: right now this fails to identify self-permutations (without any matching noise) because true self-comparisons are zero but some cross-comparisons are negative!
        # should consider pairwise consequences.
        # should consider what happens if the second choice for i_fix was chosen instead
        
        # so we fix this by editing other structures
        # first clean order matrix
        # i_fix gets its top pick, which is j_fix!
        j_fix <- rmsd_order[[ i_fix ]][1]
        # and second choice is this
        j_fix2 <- rmsd_order[[ i_fix ]][2]

        # the calculation here is to find the next best i that would have chosen j_fix
        # get the top "j" choices for every column as of now
        js <- sapply( rmsd_order, function(x) x[1] ) # get first element
        # these other "i"s could claim j_fix as their choice (includes i_fix)
        i_fix2 <- which( js == j_fix )
        # remove i_fix to see what else is there
        i_fix2 <- setdiff( i_fix2, i_fix )
        # we usually start with ties, but as we resolve them we can be left with singletons, in which case nothing else has to be done
        if ( length( i_fix2 ) > 0 ) {
            # if there's more than one choice, pick the one with the next best score
            if ( length( i_fix2 ) > 1 )
                i_fix2 <- i_fix2[ which.min( scores[ i_fix2 ] ) ]
            
            # the last thing that matters here is that is' second choice be j_fix2 as well
            # if not, then the pairwise test doesn't make sense
            if ( rmsd_order[[ i_fix2 ]][2] == j_fix2 ) {
                # now compare the combined pair scores both ways
                s1 <- rmsd_mat[ i_fix, j_fix ] + rmsd_mat[ i_fix2, j_fix2 ]
                s2 <- rmsd_mat[ i_fix, j_fix2 ] + rmsd_mat[ i_fix2, j_fix ]
                # the first score is what we get if we keep j_fix as the next one to fix
                # so we should take the other choice if s2 is smaller
                if ( s2 < s1 )
                    j_fix <- j_fix2
            }
        }
        
        # fix the best choice (after the potential pairwise analysis)
        rmsd_order[[ i_fix ]] <- j_fix
        
        # now j_fix gets removed as option from all other cases
        for ( i in 1 : k ) {
            rmsd_order_i <- rmsd_order[[ i ]]
            # skip cases that have been fixed already
            if ( length( rmsd_order_i ) == 1 )
                next
            # remove case that was just fixed elsewhere
            rmsd_order[[ i ]] <- setdiff( rmsd_order_i, j_fix )
        }
        # now update indexes_done by appending j_fix
        indexes_done <- c( indexes_done, j_fix )
    }
    
    # we have fixed all pairs now!
    # extract desired order map
    indexes <- unlist( rmsd_order ) # each has one element!
    # check for uniqueness, to validate algorithm
    stopifnot ( length( indexes ) == length( unique( indexes ) ) )
    # done!
    return( indexes )
}

# computes root mean squared error between two admixture matrix
# only validates dimensions
# doesn't require Q's to be strictly valid (i.e. non-negative and have rows sum to 1), will return RMSDs in those cases too!
rmsd_Q <- function( Q, Q2 ) {
    # validations
    if ( !is.matrix( Q ) )
        stop( '`Q` must be a matrix!' )
    if ( !is.matrix( Q2 ) )
        stop( '`Q2` must be a matrix!' )
    
    # dimensions
    n <- nrow(Q)
    k <- ncol(Q)

    if ( nrow( Q2 ) != n )
        stop( '`Q2` has ', nrow( Q2 ), ' rows, expected ', n )
    if ( ncol( Q2 ) != k )
        stop( '`Q2` has ', ncol( Q2 ), ' columns, expected ', k )

    # "norm" does all the hard work efficiently
    return( norm( Q - Q2, "F" ) / sqrt( n * k ) )
}



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
