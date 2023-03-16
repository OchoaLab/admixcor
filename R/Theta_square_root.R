# uses EVD to calculate a matrix square root
# it is truncated to the top K dimensions, so it's not an exact square root
Theta_square_root <- function( Theta, K ) {
    ev <- eigen( Theta )
    evec <- ev$vectors[ , 1L:K ]
    eval <- diag( sqrt( ev$values[ 1L:K ] ) )
    ThetaSR <- evec %*% eval
    return( ThetaSR )
}
