update_Q <- function( ThetaSR, L, R, I ) {
    
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

    # break out of here if we encounter a singular case, and re-draw everything outside
    # https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
    if ( inherits(try(solve(L), silent = TRUE), "try-error") )
        return( list( L_singular = TRUE ) )

    # calculate some matrices shared by all individuals
    # L guaranteed to be non-singular if we've gotten this far
    D <- solve( t( L ) )
    
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
        
        # perform the desired calculation, save column in desired place
        Q[ i, ] <- quadprog::solve.QP( D, d, C, c, meq, factorized = TRUE )$solution
    }

    return( list( Q = Q, L_singular = FALSE ) )
}
