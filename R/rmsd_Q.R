# computes root mean squared error between two admixture matrix
# only validates dimensions
# doesn't require Q's to be strictly valid (i.e. non-negative and have rows sum to 1), will return RMSDs in those cases too!

#' @export
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
