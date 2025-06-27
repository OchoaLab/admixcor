update_Q <- function( ThetaSR, L, R, I ) {

    # OLD derivation notes, a bit confusing given latest quadprog setup, but there is a connection
    # first set up the problem to resemble a linear regression for Q, of the form Ax=b
    # Q %*% ( L %*% R ) = ThetaSR
    # transpose to best match desired eq setup
    # t( L %*% R ) %*% t( Q ) = t( ThetaSR )
    #B <- t( ThetaSR ) # no need to explicitly transpose it
    # in this case A is the same for all rows
    #A <- t( L %*% R )
    # this has dim(x) = c(K, n)
    d <- tcrossprod( L %*% R, ThetaSR )
    
    # calculate some matrices shared by all individuals
    # L guaranteed to be non-singular if we've gotten this far (was tested outside)
    D <- solve( t( L ) )
    
    # start constructing output Q in parts
    n <- nrow( ThetaSR )
    K <- ncol( ThetaSR )

    # build constraint matrix (called Amat in solve.QP)
    # NOTE: these are the same across iterations so could be built outside loop to make more efficient!
    C <- cbind( 1, I, -I )
    # and constraint vector (called bvec in solve.QP)
    c <- c( 1, rep.int( 0, K ), rep.int( -1, K ) )
    # only first one is equality constraint
    meq <- 1

    # solve each row of Q separately
    Q <- matrix( 0, n, K )
    for ( i in 1L : n ) {
        # perform the desired calculation, save column in desired place
        # NOTE: is `d[,i]` missing a minus sign??? (that's what wikipedia suggests).  NO, quadprog has minus sign built in, but wikipedia doesn't
        Q[ i, ] <- quadprog::solve.QP( D, d[ , i ], C, c, meq, factorized = TRUE )$solution
    }

    return( Q )
}
