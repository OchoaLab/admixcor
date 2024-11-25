# L initialized to I
# R is initialized to I

initialize <- function(
                       Theta,
                       K,
                       n,
                       L_type = c('diagrandom', 'random')
                       ) {
    # process options
    L_type <- match.arg( L_type )

    # initialize Q randomly
    Q <- matrix( stats::runif( n * K ), ncol = K, nrow = n )
    Q <- Q / rowSums( Q )

    # initialize L
    if ( L_type == 'diagrandom' ) {
        L <- diag( stats::runif( K ), K, K )
    } else if ( L_type == 'random' ) {
        # random values between 0 and 1/sqrt(K), ensure Psi is between 0 and 1
        L <- matrix( stats::runif( K^2 ), K, K ) / sqrt(K)
        # ensure this is like cholesky
        L[ upper.tri( L ) ] <- 0
    } else
        stop( '`L_type` not implemented: ', L_type )
    
    # initialize R to identity matrix
    R <- diag( 1, K, K )

    # return everything!
    return(
        list(
            Q = Q,
            L = L,
            R = R
        )
    )
}

