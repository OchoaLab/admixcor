# uses EVD to calculate a matrix square root
# it is truncated to the top K dimensions, so it's not an exact square root
theta_square_root <- function( Theta, K ) {
    evd <- eigen( Theta )
    theta_square_root_from_trunc_evd( evd, K )
}

# version that assumes truncated EVD has been calculated outside
# K is for further truncation (as we might have calculated EVD on the largest K but want to consider lower K now); for now not optional
theta_square_root_from_trunc_evd <- function( evd, K ) {
    # some quick checks before the main computation
    if ( !is.list( evd ) )
        stop( '`evd` must be a list!' )
    names_missing <- setdiff( c('vectors', 'values'), names( evd ) )
    if ( length( names_missing ) > 0 )
        stop( '`evd` is missing these named elements: ', toString( names_missing ) )
    # TODO: have not checked that the data has at least K vectors/values!
    
    evec <- evd$vectors[ , 1L:K ]
    eval <- diag( sqrt( evd$values[ 1L:K ] ) )
    ThetaSR <- evec %*% eval
    return( ThetaSR )
}
