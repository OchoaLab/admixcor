#' @export
objective2 <- function( Theta, Q, Psi, alpha ) {
    # original version (not linearized)
    ## # reconstruct Psi
    ## Psi <- tcrossprod( L )
    ## # calculate main objective term
    ## O <- norm( Theta - tcrossprod( Q %*% Psi, Q ), "F" )^2

    # original objective
    O <- norm( Theta - tcrossprod( Q %*% Psi, Q ), "F" )^2

    # calculate penalty terms
    pPsi <- alpha * sum( diag( Psi ) )
    # sum terms
    O_total <- O + pPsi
    # return sum and parts
    return( c( O_total, O, pPsi ) )
}
