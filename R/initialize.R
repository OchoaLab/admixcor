initialize <- function(
                       Theta,
                       K,
                       n,
                       L_type = c('diagrandom', 'diageven', 'diagevensqrt')
                       ) {
    # process options
    L_type <- match.arg( L_type )

    # initialize Q randomly
    Q <- matrix( stats::runif( n * K ), ncol = K, nrow = n )
    Q <- Q / rowSums( Q )

    # initialize L
    if ( L_type == 'diagrandom' ) {
        L <- diag( stats::runif( K ), K, K )
    } else if ( L_type == 'diageven' ) {
        L <- diag( ( 1 : K ) / K, K, K )
    } else if ( L_type == 'diagevensqrt' ) {
        L <- diag( sqrt( ( 1 : K ) / K ), K, K )
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

