initialize <- function(
                       Theta,
                       K,
                       n
                       ) {
    # initialize Q randomly
    Q <- matrix( stats::runif( n * K ), ncol = K, nrow = n )
    Q <- Q / rowSums( Q )

    # initialize L, "diag even sqrt" ensures even spacing of Psi diagonal
    L <- diag( sqrt( ( 1 : K ) / K ), K, K )
    
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

