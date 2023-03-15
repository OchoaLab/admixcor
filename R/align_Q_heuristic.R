# UNUSED, an old failed experiment...

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
