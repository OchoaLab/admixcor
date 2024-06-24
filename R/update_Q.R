update_Q <- function( ThetaSR, L, R, delta, I, algorithm = c('original', 'nnls', 'bvls', 'glmnet', 'quadprog', 'quadprog-compact'), vertex_refine = FALSE ) {
    # process options
    algorithm <- match.arg( algorithm )

    if ( algorithm == 'original' ) {
        # Q <- tcrossprod( ThetaSR, R ) %*% MASS::ginv( L ) # no Q regularization
        Q <- tcrossprod( ThetaSR, L %*% R ) %*% MASS::ginv( tcrossprod( L ) + delta * I ) # regularized
    
        # project Q to simplex
        Q <- t( apply( Q, 1L, projsplx ) )
        ## # old projection hack, did not perform as well
        ## for (j in 1:n) {
        ##     Q[j,] <- ifelse( Q[j,] < 0.0, Q[j,] - min(Q[j,]) + 0.0001, Q[j,] )
        ## }
        ## Q <- Q / rowSums(Q)

        # this is only case not already covered by vertex_refine in another loop below
        # a copy of what we have elsewhere...
        if ( vertex_refine ) {
            A <- t( L %*% R )
            D <- tcrossprod( L ) + delta * I
            n <- nrow( ThetaSR )
            for ( i in 1L : n ) {
                b <- ThetaSR[ i, ]
                d <- drop( crossprod( A, b ) )
                data <- update_Q_vertex_scan( D, d )
                obj_x <- obj_quadprog( D, d, Q[ i, ] )
                if ( data$obj < obj_x )
                    Q[ i, ] <- data$q
            }
        }
    } else if ( algorithm == 'nnls' || algorithm == 'bvls' || algorithm == 'glmnet' || algorithm == 'quadprog' || algorithm == 'quadprog-compact' ) {
        # uses non-negative (nnls) or bounded variable least squares (bvls), which are close to our problem without regularization; glmnet adds regularization; all are improvement over original which has no constraints at all (at first), but none of them (except quadprog) directly restrict to simplex so projection is still necessary afterwards
        
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
        if ( algorithm == 'quadprog' || algorithm == 'quadprog-compact' || vertex_refine ) {
            # calculate some matrices shared by all individuals
            # NOTE: delta shoulds be 1/2 I think!!! (to agree with our standard objective)
            #D <- crossprod( A ) + delta * I # equivalent to below
            D <- tcrossprod( L ) + delta * I
            # build constraint matrix (called Amat in solve.QP)
            # NOTE: these are the same across iterations so could be built outside loop to make more efficient!
            if ( algorithm == 'quadprog' ) {
                C <- cbind( 1, I, -I )
            } else if ( algorithm == 'quadprog-compact' ) {
                # compact version has weird notation, but might perform better given the sparsity of our problem
                Cind <- matrix( 0L, nrow = K + 1L, ncol = 2L * K + 1L )
                Cmat <- matrix( 0L, nrow = K, ncol = 2L * K + 1L )
                # add first constraint, only dense one
                Cind[ 1L, ] <- 1L # first value is number of non-zero values, one for all but the first case (to overwrite next)
                Cind[ 1L, 1L ] <- K # first value is number of non-zero values
                Cind[ 2L : ( K + 1L ), 1L ] <- 1L : K # list all indexes for frows contiguously
                Cmat[ , 1L ] <- 1L # coefficents
                # now list every sparse constraint, each block as needed
                Cind[ 2L, 2L : (K + 1L) ] <- 1L : K
                Cind[ 2L, (K + 2L) : (2L * K + 1L) ] <- 1L : K
                # all coefficients are 1 for first block
                Cmat[ 1L, 2L : (K + 1L) ] <- 1L
                # and -1 for last block
                Cmat[ 1L, (K + 2L) : (2L * K + 1L) ] <- -1L
            }
            # and constraint vector (called bvec in solve.QP)
            c <- c( 1, rep.int( 0, K ), rep.int( -1, K ) )
            # only first one is equality constraint
            meq <- 1
        }
        # solve each row of Q (i.e. column of t(Q)) separately
        for ( i in 1L : n ) {
            # get corresponding row of ThetaSR (column of its transpose)
            b <- ThetaSR[ i, ]
            if ( algorithm == 'quadprog' || algorithm == 'quadprog-compact' || vertex_refine ) {
                # unfortunately this varies per individual because b varies
                # NOTE: is it missing a minus sign??? (that's what wikipedia suggests, have to check against quadprog package).  NO, quadprog has minus sign built in, but wikipedia doesn't
                d <- drop( crossprod( A, b ) )
                # NOTE: consider compact (sparse) version of problem!  Also pre-factorizing D into "R-inv", to share for all individuals (but it will change across iterations)
                if ( algorithm == 'quadprog' ) {
                    x <- quadprog::solve.QP( D, d, C, c, meq )$solution
                } else if ( algorithm == 'quadprog-compact' ) {
                    x <- quadprog::solve.QP.compact( D, d, Cmat, Cind, c, meq )$solution
                }
            }
            if ( algorithm %in% c('nnls', 'bvls', 'glmnet') ) {
                # solve system of equations with non-negative constraint, or both [0,1] bounds
                if ( algorithm == 'nnls' ) {
                    x <- nnls::nnls( A, b )$x
                } else if ( algorithm == 'bvls' ) {
                    x <- bvls::bvls( A, b, rep.int( 0, K ), rep.int( 1, K ) )$x
                } else if ( algorithm == 'glmnet' ) {
                    # `alpha = 0`: use "ridge" penalty
                    # reduced default tolerance of 1e-7 because exact solutions were not good enough
                    x <- glmnet::glmnet( A, b, alpha = 0, lambda = delta, intercept = FALSE, lower.limits = 0, upper.limits = 1, thresh = 1e-20 )$beta[,1]
                }
                # now project to simplex (fixes upper bound cap for nnls, only case that doesn't enforce it directly)
                x <- projsplx( x )
            }

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

    }

    return( Q )
}
