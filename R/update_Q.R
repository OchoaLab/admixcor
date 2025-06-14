update_Q <- function( ThetaSR, L, R, delta, I ) {
    
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

    # troubleshooting stats
    # remember if L was singular
    # https://stackoverflow.com/questions/24961983/how-to-check-if-a-matrix-has-an-inverse-in-the-r-language
    L_singular <- inherits(try(solve(L), silent = TRUE), "try-error")

    # break out of here if we encounter a singular case, and re-draw everything outside
    if ( L_singular )
        return( list( L_singular = L_singular ) )

    # calculate some matrices shared by all individuals
    #D <- crossprod( A ) + delta * I # equivalent to below when delta!=0
    # if delta=0, construct a "factored version" that's supposed to be faster
    factorized <- FALSE
    if ( delta == 0 ) {
        # L guaranteed to be non-singular if we've gotten this far
        D <- solve( t( L ) )
        factorized <- TRUE
    } else {
        D <- tcrossprod( L ) + delta * I
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
        
        # perform the desired calculation, save column in desired place
        Q[ i, ] <- quadprog::solve.QP( D, d, C, c, meq, factorized = factorized )$solution
    }

    return( list( Q = Q, L_singular = L_singular ) )
}
