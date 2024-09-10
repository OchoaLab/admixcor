# update L (square root of Psi)
update_L <- function( ThetaSR, Q, R, gamma = 0, algorithm = c('glmnet', 'bvls'), fix_L = FALSE ) {
    # process options
    algorithm <- match.arg( algorithm )

    # uses non-negative (nnls) or bounded variable least squares (bvls), which are equivalent to our problem without regularization; glmnet in theory solves our problem with regularization!
    # nnls solves exactly the requirement that the lower bound be zero, and exactly forces lower triangle to be zero during optimization, but does not set upper bound directly, we set strict cap of 1 as before

    # first set up the problem to resemble a linear regression for L, of the form Ax=b
    
    # start constructing the output L
    K <- ncol( Q )
    L <- matrix( 0, K, K )
    
    # https://en.wikipedia.org/wiki/Kronecker_product#Matrix_equations
    # first version assumes L is a full matrix, and the notation is quite standard for Kronecker product (%x%) vec trick (as.vector(matrix))
    A <- t( R ) %x% Q
    b <- as.vector( ThetaSR )
    # now, half of L is zero, let's make sure that's reflected in the problem statement (so that constraint is enforced)
    indexes <- as.vector( if ( fix_L ) lower.tri( L, diag = TRUE ) else upper.tri( L, diag = TRUE ) )
    A <- A[ , indexes ]
    # now get solutions!
    # solve system of equations with non-negative constraint, or both [0,1] bounds
    if ( algorithm == 'bvls' ) {
        if ( gamma != 0 )
            stop( '`L_algorithm == "bvls"` requires `gamma = 0`' )
        # length of x, needed to set bounds
        lx <- K * ( K + 1 ) / 2
        x <- bvls::bvls( A, b, rep.int( 0, lx ), rep.int( 1, lx ) )$x
    } else if ( algorithm == 'glmnet' ) {
        if ( K == 1L ) {
            # glmnet doesn't work to solve a single variable, super annoying!
            # I wrote my own solver for that nearly trivial case
            x <- glmnet_one_bounded( A, b, gamma )
        } else {
            # `alpha = 0`: use "ridge" penalty
            # reduced default tolerance of 1e-7 because exact solutions were not good enough
            x <- glmnet::glmnet( A, b, alpha = 0, lambda = gamma, intercept = FALSE, lower.limits = 0, upper.limits = 1, thresh = 1e-20 )$beta[,1]
        }
    }
    # lastly, reform x into L, return that!
    # (L was already a zero matrix, so the rest is all set!)
    L[ indexes ] <- x
    ## if ( fix_L ) {
    ##     L[ lower.tri( L, diag = TRUE ) ] <- x
    ## } else {
    ##     L[ upper.tri( L, diag = TRUE ) ] <- x
    ## }
    
    return( L )
}

glmnet_one_bounded <- function( a, b, gamma ) {
    # this is very silly, but glmnet won't work if there's a single variable to solve for
    # so let's do it explicitly
    # this is objective, for reference:
    # O = ( a * x - b ) %*% t( a * x - b ) + gamma * x^2
    # this is unbounded solution
    #x <- ( b %*% a ) / ( crossprod( a ) + gamma )
    x <- crossprod( a, b ) / ( crossprod( a ) + gamma )

    # if solution is in bounds, that's the desired solution!
    if ( 0 < x && x < 1 )
        return( x )
    
    # else, if solution is out of bounds, for this single variable case only other solutions are edge cases, i.e. 0 and 1
    obj0 <- crossprod( b )
    obj1 <- crossprod( a - b ) + gamma
    return( if ( obj0 > obj1 ) 0 else 1 )
}
