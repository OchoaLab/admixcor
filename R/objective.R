#' @export
objective <- function( ThetaSR, Q, L, R, gamma ) {
    ## # original version (not linearized)
    ## # reconstruct Psi
    ## Psi <- tcrossprod( L )
    ## # calculate main objective term
    ## O <- norm( Theta - tcrossprod( Q %*% Psi, Q ), "F" )^2

    # calculate linearized objective, the one we're actually optimizing!
    O <- norm( ThetaSR - Q %*% L %*% R, "F" )^2

    # calculate penalty terms
    pL <- gamma * norm( L, "F" )^2
    # sum terms
    O_total <- O + pL
    # return sum and parts
    return( c( O_total, O, pL ) )
}
