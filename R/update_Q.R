update_Q <- function( ThetaSR, L, R, delta, I, vertex_refine = FALSE ) {

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
    #D <- crossprod( A ) + delta * I # equivalent to below
    D <- tcrossprod( L ) + delta * I
    # hack-correct D if needed, a minimal alteration so we don't have to repeat runs that previously succeeded
    if ( delta == 0 ) {
        eps <- 1e-2 # sqrt( .Machine$double.eps )
        indexes <- diag( D ) < eps
        if ( any( indexes ) ) {
            # set things that are too small to the smallest non-zero value that is guaranteed to pass (make D posdef)
            diag( D )[ indexes ] <- eps
        }
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
        # NOTE: consider compact (sparse) version of problem!  Also pre-factorizing D into "R-inv", to share for all individuals (but it will change across iterations)
        x <- quadprog::solve.QP( D, d, C, c, meq )$solution
        
        if ( vertex_refine ) {
            # test vertices (all K unadmixed cases), see if the best one outperform the solution of the regular algorithms (if so, replace solution)
            data <- update_Q_vertex_scan( D, d )
            # now calculate objective of current solution
            obj_x <- obj_quadprog( D, d, x )
            # replace if vertex objective is better
            if ( data$obj < obj_x )
                x <- data$q
        }
        # save column in desired place
        Q[ i, ] <- x
    }

    return( Q )
}
