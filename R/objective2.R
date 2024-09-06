#' @export
objective2 <- function( Theta, Q, Psi, alpha, beta ) {
    ## # original version (not linearized)
    ## # reconstruct Psi
    ## Psi <- tcrossprod( L )
    ## # calculate main objective term
    ## O <- norm( Theta - tcrossprod( Q %*% Psi, Q ), "F" )^2

    # calculate linearized objective, the one we're actually optimizing!
    O <- norm( Theta - tcrossprod( Q %*% Psi, Q ), "F" )^2

    # calculate penalty terms
    pPsi <- alpha * sum( diag( Psi ) )
    pQ <- delta * norm( Q, "F" )^2
    # sum terms
    O_total <- O + pPsi + pQ
    # return sum and parts
    return( c( O_total, O, pPsi, pQ ) )
}
