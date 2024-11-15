# `tol` is lowest value allowed (should be negative to have tolerance, the intolerant version is zero)
#' @export
stretch_Q <- function( Q, shrink = TRUE, ties_none = FALSE, tol = -0.01, n_iter = -log2( .Machine$double.eps ) ) {
    # need the identity matrix here
    I <- diag( ncol( Q ) )
    
    # find the individual with the maximum of each ancestry
    # this special function identifies ties and handles them, important to avoid singularity further below
    indexes <- vertex_inds( Q, ties_none = ties_none )
    # indexes can have missing values when no suitable vertex was identified, leave those cases unstretched
    S_inv <- I # to allow for unstretched cases by default
    indexes2 <- !is.na( indexes )
    # this is the matrix of current vertices Qv to use to stretch, which is also the inverse of the stretch matrix
    S_inv[ indexes2, ] <- Q[ indexes[ indexes2 ], ]
    # now get its inverse, which is the stretching transformation
    # in very malformed data this is not invertible despite our best attempts (in one example an entire ancestry was just not predicted, that column was all zeroes), a generalized inverse is guaranteed to work!
    S <- MASS::ginv( S_inv )
    # apply to Q, let's see if it stayed in range or not
    Q2 <- Q %*% S
    # the implicit alpha, to return if things are good or to modify otherwise
    alpha <- 1
    
    if ( min( Q2 ) < tol && shrink ) {
        # mark where we were with current alpha
        alpha_min <- 0
        alpha_max <- 1
        # first time it's too high
        high <- TRUE
        # perform a continuous binary search
        # the number of iterations is fixed, set by machine precision
        for ( i in 1 : n_iter ) {
            if ( high ) {
                # if we're too high, current alpha becomes the max to consider from here on
                alpha_max <- alpha
                # for next alpha, consider midpoint between current alpha and minimum
                alpha <- ( alpha + alpha_min ) / 2
            } else {
                # if this was too low, then current alpha becomes the min to consider from here on
                alpha_min <- alpha
                # for next alpha, consider midpoint between current alpha and maximum
                alpha <- ( alpha + alpha_max ) / 2
            }
            # try the new alpha
            Sa <- alpha * S + ( 1 - alpha ) * I
            Q2 <- Q %*% Sa
            # decide on direction to take next
            high <- min( Q2 ) < tol 
        }
        # after we're done with this loop, replace S with the best Sa
        S <- Sa
        # S_inv has to be recalculated now
        S_inv <- MASS::ginv( S )
        # Q2 and alpha are already as desired
    }
    
    # return the updated Q, and also the transformation matrices for external tests
    return( list (
        Q = Q2,
        S = S,
        S_inv = S_inv,
        alpha = alpha
    ) )
}

