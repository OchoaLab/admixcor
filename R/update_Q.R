update_Q <- function( ThetaSR, L, R, delta, I, delta2 = 1e-6 ) {
    
    # first set up the problem to resemble a linear regression for Q, of the form Ax=b
    # Q %*% ( L %*% R ) = ThetaSR
    # transpose to best match desired eq setup
    # t( L %*% R ) %*% t( Q ) = t( ThetaSR )
    #B <- t( ThetaSR ) # no need to explicitly transpose it
    # in this case A is the same for all rows
    A <- t( L %*% R )
    # start constructing output Q in parts
    n <- nrow( ThetaSR )
    K <- ncol( ThetaSR )
    Q <- matrix( 0, n, K )

    # calculate some matrices shared by all individuals
    #D <- crossprod( A ) + delta * I # equivalent to below when delta!=0
    # if delta=0, construct a "factored version" that solves some singularity problems only observed in that case
    # in fact, L is also singular, so try a generalized inverse in those cases?
    factorized <- FALSE
    if ( delta == 0 ) {
        D <- MASS::ginv( t( L ) )
        factorized <- TRUE
        # calculate a mildly regularized version in case of errors, which sadly occur somewhat frequently
        D2 <- tcrossprod( L ) + delta2 * I
    } else {
        D <- tcrossprod( L ) + delta * I
        # errors never happen in this case, but this has to be defined
        D2 <- D
    }
    
    # build constraint matrix (called Amat in solve.QP)
    # NOTE: these are the same across iterations so could be built outside loop to make more efficient!
    C <- cbind( 1, I, -I )
    # and constraint vector (called bvec in solve.QP)
    c <- c( 1, rep.int( 0, K ), rep.int( -1, K ) )
    # only first one is equality constraint
    meq <- 1

    # TODO: below each little `d` could probably be precomputed as a matrix, but we'll save testing that idea for later
    
    # solve each row of Q (i.e. column of t(Q)) separately
    for ( i in 1L : n ) {
        # get corresponding row of ThetaSR (column of its transpose)
        b <- ThetaSR[ i, ]

        # unfortunately this varies per individual because b varies
        # NOTE: is it missing a minus sign??? (that's what wikipedia suggests, have to check against quadprog package).  NO, quadprog has minus sign built in, but wikipedia doesn't
        d <- drop( crossprod( A, b ) )
        
        # try the problem we want to solve, though when delta==0 we may have unexpected errors
        x <- try(
            quadprog::solve.QP( D, d, C, c, meq, factorized = factorized )$solution,
            silent = TRUE
        )
        # if there was an error, use this mildly regularized version never causes errors!
        if ( inherits( x, "try-error" ) )
            x <- quadprog::solve.QP( D2, d, C, c, meq )$solution
        # save column in desired place
        Q[ i, ] <- x
    }

    return( Q )
}
