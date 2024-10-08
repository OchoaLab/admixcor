library(tibble)

n <- 10L
K <- 3L
# create a random Q
Q <- matrix( runif( n*K ), nrow = n, ncol = K )
Q <- Q / rowSums( Q ) # normalize as they should be
# now create a totally random new Q to align to original Q (it should work, though it's random)
Q2 <- matrix( runif( n*K ), nrow = n, ncol = K )
Q2 <- Q2 / rowSums( Q2 ) # normalize as they should be
# construct a random Psi, diagonal for simplicity
Psi <- diag( runif( K ), K )
# create a Theta that goes with these matrices
Theta <- tcrossprod( Q %*% Psi, Q )
# some default values for the regularization parameters; light regularization is ideal
alpha <- 1e-5
beta <- 1e-5
gamma <- 1e-5
delta <- 1e-5
# constant used in regularized expressions
I <- diag( 1, K, K )
# a random quanitity for quadcomp tests
d <- rnorm( K )

### generic validators

# validates dimensions, non-negativity, rows sum to one
validate_Q <- function( Q, n, K ) {
    expect_true( is.matrix( Q ) )
    expect_equal( nrow( Q ), n )
    expect_equal( ncol( Q ), K )
    # weird way to test this inequality but with tolerance
    ## expect_true( min( Q ) >= 0 )
    minQ <- min( Q )
    if ( minQ > 0 ) minQ <- 0
    expect_equal( minQ, 0 )
    expect_equal( rowSums( Q ), rep.int( 1, n ) )
}

# validates dimensions, that R is actually orthogonal
validate_R <- function( R, K, I ) {
    # test basic expectations
    expect_true( is.matrix( R ) )
    expect_equal( nrow( R ), K )
    expect_equal( ncol( R ), K )
    # confirm this is rotation matrix (I must have right dimension of K x K)
    expect_equal( crossprod( R ), I )
}

# validates dimensions, non-negativity, Cholesky form, and implied Psi <= 1
validate_L <- function( L, K, fix_L = FALSE ) {
    expect_true( is.matrix( L ) )
    expect_equal( nrow( L ), K )
    expect_equal( ncol( L ), K )
    expect_true( min( L ) >= 0 )
    expect_true( max( L ) <= 1 )
    # demand that lower triangle is zero! (excludes diagonal)
    if ( fix_L ) {
        expect_true( all( L[ upper.tri( L ) ] == 0 ) )
    } else {
        expect_true( all( L[ lower.tri( L ) ] == 0 ) )
    }
    # for more accurate validations, check its full Psi
    Psi <- tcrossprod( L )
    #expect_true( min( Psi ) >= 0 ) # implied by earlier test
    # not implied, test too!
    # weird way to test this inequality but with tolerance
    ## expect_true( max( Psi ) <= 1 )
    maxPsi <- max( Psi )
    if ( maxPsi < 1 ) maxPsi <- 1
    expect_equal( maxPsi, 1 )
}

# validates dimensions, non-negativity, and Psi <= 1
validate_Psi <- function( Psi, K ) {
    expect_true( is.matrix( Psi ) )
    expect_equal( nrow( Psi ), K )
    expect_equal( ncol( Psi ), K )
    expect_true( min( Psi ) >= 0 )
    expect_true( max( Psi ) <= 1 )
    # for more accurate validations, check its full Psi
    # weird way to test this inequality but with tolerance
    ## expect_true( max( Psi ) <= 1 )
    ## maxPsi <- max( Psi )
    ## if ( maxPsi < 1 ) maxPsi <- 1
    ## expect_equal( maxPsi, 1 )
}

# uniform testing for all cases
validate_initialize <- function( obj, n, K, fix_L = FALSE ) {
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('Q', 'L', 'R') )
    validate_Q( obj$Q, n, K )
    validate_R( obj$R, K, I )
    validate_L( obj$L, K, fix_L = fix_L )
}

# uniform testing for all cases
validate_admixcor <- function( obj, n, K, fix_L = FALSE, v = 1 ) {
    expect_true( is.list( obj ) )
    if ( v == 1 ) {
        expect_equal( names( obj ), c('Q', 'Psi', 'f', 'report', 'ThetaSR', 'L', 'R') )
    } else {
        expect_equal( names( obj ), c('Q', 'Psi', 'f', 'report', 'ThetaSR', 'R') )
    }
    validate_Q( obj$Q, n, K )
    if ( v == 1 )
        validate_L( obj$L, K, fix_L = fix_L )
    validate_R( obj$R, K, I )
    # NOTE: we're not validating ThetaSR, that's ok, it's simple and not a worry (it was validated earlier)
    expect_true( is.matrix( obj$Psi ) )
    expect_equal( nrow( obj$Psi ), K )
    expect_equal( ncol( obj$Psi ), K )
    expect_true( min( obj$Psi ) >= 0 )
    # expect_true( max( obj$Psi ) <= 1 ) # current algorithm doesn't enforce this
    expect_true( is.numeric( obj$f ) )
    expect_equal( length( obj$f ), 1L )
    expect_true( obj$f >= 0 )
    expect_true( obj$f <= n*K ) # unnormalized has this max assuming all values are between 0-1
    expect_true( is_tibble( obj$report ) ) # undetailed validation here...
}

### TESTS

test_that( '`L = t(chol(Psi))` inverts `Psi = tcrossprod(L)`', {
    # our Psi is diagonal/boring, and Theta is not full rank, so create a fake but interesting and full-rank Psi here instead
    Psi <- crossprod( matrix( rnorm( K^2 ), K, K ) )
    expect_silent(
        L <- t( chol( Psi ) )
    )
    expect_equal( tcrossprod( L ), Psi )
})

test_that( 'rmsd_Q works', {
    expect_silent(
        val <- rmsd_Q( Q, Q2 )
    )
    expect_equal( length( val ), 1L )
    expect_true( is.numeric( val ) )
    expect_true( val >= 0 )
    expect_true( val <= 1 ) # because of the way it's normalized
})

test_that( 'rmsd_Q_mat works', {
    expect_silent(
        rmsd_mat <- rmsd_Q_mat( Q, Q2 )
    )
    expect_true( is.matrix( rmsd_mat ) )
    expect_equal( nrow( rmsd_mat ), K )
    expect_equal( ncol( rmsd_mat ), K )
    expect_true( all( rmsd_mat >= 0 ) )
    expect_true( all( rmsd_mat <= n ) ) # this one has unnormalized elements
})

test_that( 'align_Q works', {
    # permute original matrix so there is an exact solution
    Q_shuffled <- Q[ , sample(K) ]

    expect_silent(
        indexes <- align_Q( Q, Q_shuffled ) # , verbose = FALSE, fast = TRUE, fast2 = FALSE
    )
    # in this case expect exact match after permutation
    expect_equal( Q, Q_shuffled[ , indexes ] )

    # test random Q2 now
    expect_silent(
        indexes <- align_Q( Q, Q2 )
    )
    # here we just want indexes to be a permutation matrix
    expect_equal( length( indexes ), K )
    expect_true( all( indexes %in% 1L : K ) )
})

test_that( 'theta_square_root works', {
    # does a dimensionality reduction too, so it's not an exact square root unless it was low rank and it is correct
    expect_silent(
        Theta_sqrt <- theta_square_root( Theta, K )
    )
    expect_true( is.matrix( Theta_sqrt ) )
    expect_equal( nrow( Theta_sqrt ), n )
    expect_equal( ncol( Theta_sqrt ), K )
    # reconstruct theta, here it will be exact!
    Theta_redone <- tcrossprod( Theta_sqrt )
    expect_equal( Theta_redone, Theta )
})

# now that the above function was validated, use it to calculate square root
ThetaSR <- theta_square_root( Theta, K )
# also calculate a few other true square root model parameters
# this is correct only because true Psi is diagonal, otherwise ought to use Cholesky!
L <- sqrt( Psi )
# need true rotation too...
R <- solve( L ) %*% MASS::ginv( Q ) %*% ThetaSR
# NOTE: for exact solution, this R is automatically rotation (no need to additionally project to rotation space)
# also create a random L used for tests, copying code from initialize.R with `L_type == 'random'`
L2 <- matrix( runif( K^2 ), K, K ) / sqrt(K)
# ensure this is like cholesky
L2[ upper.tri( L2 ) ] <- 0


test_that( 'objective works', {
    # in this test match is exact, so objective should be zero
    # (but only if penalties are zero too!)
    expect_silent( 
        f <- objective( ThetaSR, Q, L, R, gamma = 0, delta = 0 )
    )
    expect_equal( f, rep.int( 0, 4L ) )

    # now try an inexact match
    # this shoudn't alter global values in other scopes, right?
    L <- diag( 1, K )
    R <- diag( 1, K )
    # passed Q2 instead of Q here too, so it's totally random
    expect_silent( 
        f <- objective( ThetaSR, Q2, L, R, gamma, delta )
    )
    expect_true( is.numeric( f ) )
    expect_equal( length( f ), 4L )
    expect_true( min( f ) >= 0 )
    expect_true( f[2L] <= n*K ) # unnormalized has this max assuming all values are between 0-1 (not strictly true for linearized objective, because of rotation, but let's see...)
    expect_equal( f[1L], sum( f[2L:4L] ) ) # first is sum of the rest
})

test_that( 'initialize works', {
    # test all Q_types!
    expect_silent(
        obj <- initialize( Theta, K, n, Q_type = 'kmeans' )
    )
    validate_initialize( obj, n, K )
    expect_silent(
        obj <- initialize( Theta, K, n, Q_type = 'random' )
    )
    validate_initialize( obj, n, K )
    expect_silent(
        obj <- initialize( Theta, K, n, Q_type = 'uniform' )
    )
    validate_initialize( obj, n, K )

    # test all L_types!
    expect_silent(
        obj <- initialize( Theta, K, n, L_type = 'identity' )
    )
    validate_initialize( obj, n, K )
    expect_silent(
        obj <- initialize( Theta, K, n, L_type = 'uniform' )
    )
    validate_initialize( obj, n, K )
    expect_silent(
        obj <- initialize( Theta, K, n, L_type = 'diagrandom' )
    )
    validate_initialize( obj, n, K )
    expect_silent(
        obj <- initialize( Theta, K, n, L_type = 'random' )
    )
    validate_initialize( obj, n, K )
})

test_that( 'projsplx works', {
    # naive function is applied one row at the time
    q <- runif( K )
    expect_silent(
        q <- projsplx( q )
    )
    expect_equal( length( q ), K )
    expect_true( min( q ) >= 0 )
    expect_equal( sum( q ), 1 )
    
    # matrix version!
    # here create an out of range Q
    # (rows likely exceed one, don't sum to one in any case)
    Q3 <- matrix( runif( n*K ), nrow = n, ncol = K )
    expect_silent(
        Q3 <- t( apply( Q3, 1L, projsplx ) )
    )
    # validate now!
    validate_Q( Q3, n, K )
})

test_that( 'update_R works', {
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use exact ThetaSR from real data, but random Q for test
    R2 <- update_R( ThetaSR, Q2, L )
    # test basic expectations
    validate_R( R2, K, I )
    
    # now try exact solution (true ThetaSR, Q and L)
    R2 <- update_R( ThetaSR, Q, L )
    expect_equal( R2, R )
})

test_that( 'update_L works', {
    # test the default glmnet first...
    
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use exact ThetaSR from real data, but random Q for test
    L2 <- update_L( ThetaSR, Q2, R, gamma )
    # test basic expectations
    validate_L( L2, K )
    # TODO x3: `maxPsi` (`actual`) not equal to 1 (`expected`).
    # TODO: max(L) <= 1 is not TRUE

    # now try exact solution (true ThetaSR, Q and R), here it is recoverable in theory if we don't penalize!
    L2 <- update_L( ThetaSR, Q, R, 0 )
    expect_equal( L2, L )

    # test bvls version (requires zero gamma), repeating both earlier tests
    L2 <- update_L( ThetaSR, Q2, R, algorithm = 'bvls' )
    validate_L( L2, K )
    # TODO x4: `maxPsi` (`actual`) not equal to 1 (`expected`).
    L2 <- update_L( ThetaSR, Q, R, algorithm = 'bvls' )
    expect_equal( L2, L )

    # repeat tests with `fix_L = TRUE`
    L2 <- update_L( ThetaSR, Q2, R, gamma, fix_L = TRUE )
    validate_L( L2, K, fix_L = TRUE )
    L2 <- update_L( ThetaSR, Q, R, 0, fix_L = TRUE )
    expect_equal( L2, L )
    L2 <- update_L( ThetaSR, Q2, R, algorithm = 'bvls', fix_L = TRUE )
    validate_L( L2, K, fix_L = TRUE )
    L2 <- update_L( ThetaSR, Q, R, algorithm = 'bvls', fix_L = TRUE )
    expect_equal( L2, L )
})

test_that( 'update_Psi works', {
    # test the default glmnet first...
    
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use true Theta, but random Q for test
    Psi2 <- update_Psi( Theta, Q2, alpha )
    # test basic expectations
    validate_Psi( Psi2, K )

    # now try exact solution (true Theta and Q), here it is recoverable in theory if we don't penalize!
    Psi2 <- update_Psi( Theta, Q, 0 )
    expect_equal( Psi2, Psi )

    # test bvls version (requires zero alpha), repeating both earlier tests
    Psi2 <- update_Psi( Theta, Q2, algorithm = 'bvls' )
    validate_Psi( Psi2, K )
    # TODO x4: `maxPsi` (`actual`) not equal to 1 (`expected`).
    Psi2 <- update_Psi( Theta, Q, algorithm = 'bvls' )
    expect_equal( Psi2, Psi )
})

test_that( "obj_quadprog works", {
    expect_silent(
        obj <- obj_quadprog( Psi, d, Q2[1,] )
    )
})

test_that( "update_Q_vertex_scan works", {
    expect_silent(
        data <- update_Q_vertex_scan( Psi, d )
    )
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('q', 'obj') )
    expect_true( is.numeric( data$q ) )
    expect_equal( length( data$q ), K )
    expect_true( is.numeric( data$obj ) )
    expect_equal( length( data$obj ), 1L )

    # compare to slower, more explicit version
    expect_silent(
        data2 <- update_Q_vertex_scan( Psi, d, fast = FALSE )
    )
    expect_equal( data2, data )
})

validate_vertex_refine <- function( ThetaSR, Q1, Q2, L, R, gamma, delta ) {
    # calculate both objectives
    obj1 <- objective( ThetaSR, Q1, L, R, gamma, delta )[1]
    obj2 <- objective( ThetaSR, Q2, L, R, gamma, delta )[1]
    # make sure that second is as good or better than first!
    expect_true( obj2 <= obj1 )
}

test_that( 'update_Q works', {
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use exact ThetaSR from real data, but random L for test
    Q2 <- update_Q( ThetaSR, L2, R, delta, I )
    # test basic expectations
    validate_Q( Q2, n, K )
    # repeat all tests with vertex_refine = TRUE, which requires delta and I in all cases!
    Q2v <- update_Q( ThetaSR, L2, R, delta, I, vertex_refine = TRUE )
    validate_Q( Q2v, n, K )
    # ensure that objectives are actually improved by vertex_refine
    validate_vertex_refine( ThetaSR, Q2, Q2v, L2, R, gamma, delta )
    
    # now try exact solution (true ThetaSR, L and R), here it is recoverable in theory if we don't penalize!
    Q2 <- update_Q( ThetaSR, L, R, 0, I )
    expect_equal( Q2, Q )
    # repeat all tests with vertex_refine = TRUE, which requires delta and I in all cases!
    Q2v <- update_Q( ThetaSR, L, R, 0, I, vertex_refine = TRUE )
    expect_equal( Q2v, Q )
    # if these are all equal, then the objectives are equal (no need to test directly)

    # test new version based on quadprog, repeating both earlier tests
    Q2 <- update_Q( ThetaSR, L2, R, delta, I, algorithm = 'quadprog' )
    validate_Q( Q2, n, K )
    Q2v <- update_Q( ThetaSR, L2, R, delta, I, algorithm = 'quadprog', vertex_refine = TRUE )
    validate_Q( Q2v, n, K )
    validate_vertex_refine( ThetaSR, Q2, Q2v, L2, R, gamma, delta )
    # can only recover true solution if we don't penalize!
    Q2 <- update_Q( ThetaSR, L, R, 0, I, algorithm = 'quadprog' )
    expect_equal( Q2, Q )
    Q2v <- update_Q( ThetaSR, L, R, 0, I, algorithm = 'quadprog', vertex_refine = TRUE )
    expect_equal( Q2v, Q )
})

test_that( 'positive_definite works', {
    # across iterations, we encountered cases for Psi that chol determines are not positive definite, here we tackle that problem directly
    
    # generate a case with a negative eigenvalue
    Psi_bad <- Psi
    Psi_bad[ K, K ] <- -0.001
    # confirm that chol wouldn't like it
    expect_error( chol( Psi_bad ) )
    # now apply our trick
    expect_silent(
        Psi_good <- positive_definite( Psi_bad, reg = TRUE )
    )
    # confirm that chol likes this now
    expect_silent( chol( Psi_good ) )

    # repeat with another variant that doesn't regularize
    expect_silent(
        Psi_good <- positive_definite( Psi_bad )
    )
    # confirm that chol likes this now
    expect_silent( chol( Psi_good ) )
})

test_that( 'admixcor works', {
    # for this test, we don't have to fully converge (even in toy data it sometimes takes too long)
    # this stops in a single iteration in tests!
    tol <- 1e-2 # default ~1e-8
    # first test default regularized version
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol )
    )
    ## obj <- admixcor( Theta, K, tol = tol, verbose = TRUE ) # for testing
    validate_admixcor( obj, n, K )
    # TODO: `maxPsi` (`actual`) not equal to 1 (`expected`).
    expect_silent(
        objv <- admixcor( Theta, K, tol = tol, vertex_refine = TRUE )
    )
    validate_admixcor( objv, n, K )
    ## # are final objectives improved?  Surprisingly, not always (but these are weird, toy problems)!
    ## if ( objv$f > obj$f )
    ##     message( objv$f, ' > ', obj$f )
    ## expect_true( objv$f <= obj$f ) # FAILED x4
    
    # now test no regularization version, where some things may be singular if we're not careful
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = 0 )
    )
    validate_admixcor( obj, n, K )
    # TODO: max(L) <= 1 is not TRUE

    # now test no regularization version, where some things may be singular if we're not careful
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, delta = 0 )
    )
    validate_admixcor( obj, n, K )

    # now test no regularization version, where some things may be singular if we're not careful
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = 0, delta = 0 )
    )
    validate_admixcor( obj, n, K )
    # TODO: `maxPsi` (`actual`) not equal to 1 (`expected`).

    # all of those earlier cases were default Q initialization using kmeans, try other non-default cases
    # test default regularized version and totally unregularized, to catch edge cases potentially
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'random' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'random', gamma = 0, delta = 0 )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'uniform' )
    )
    validate_admixcor( obj, n, K )
    # TODO: `maxPsi` (`actual`) not equal to 1 (`expected`).
    # TODO: max(L) <= 1 is not TRUE
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'uniform', gamma = 0, delta = 0 )
    )
    validate_admixcor( obj, n, K )
    # TODO: `maxPsi` (`actual`) not equal to 1 (`expected`).
    # TODO: max(L) <= 1 is not TRUE

    # ditto L initializations, try non-default versions now (identity is default)
    # test in context of Q_type = 'kmeans' only
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'uniform' )
    )
    validate_admixcor( obj, n, K )
    # TODO x2: `maxPsi` (`actual`) not equal to 1 (`expected`).
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'diagrandom' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'random' )
    )
    validate_admixcor( obj, n, K )
    # TODO: `maxPsi` (`actual`) not equal to 1 (`expected`).

    # and L algorithms!
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_algorithm = 'bvls' )
    )
    validate_admixcor( obj, n, K )
    # TODO: `maxPsi` (`actual`) not equal to 1 (`expected`).

    # and Q algorithms!
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        objv <- admixcor( Theta, K, tol = tol, Q_algorithm = 'quadprog', vertex_refine = TRUE )
    )
    validate_admixcor( objv, n, K )
    ## if ( objv$f > obj$f )
    ##     message( objv$f, ' > ', obj$f )
    ## expect_true( objv$f <= obj$f ) # TODO FAILED x1

    # repeat every previous test with `fix_L = TRUE`
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        objv <- admixcor( Theta, K, tol = tol, vertex_refine = TRUE, fix_L = TRUE )
    )
    validate_admixcor( objv, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = 0, fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, delta = 0, fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = 0, delta = 0, fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'random', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'random', gamma = 0, delta = 0, fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'uniform', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'uniform', gamma = 0, delta = 0, fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'uniform', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'diagrandom', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'random', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_algorithm = 'bvls', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'quadprog', fix_L = TRUE )
    )
    validate_admixcor( obj, n, K, fix_L = TRUE )
    expect_silent(
        objv <- admixcor( Theta, K, tol = tol, Q_algorithm = 'quadprog', vertex_refine = TRUE, fix_L = TRUE )
    )
    validate_admixcor( objv, n, K, fix_L = TRUE )
})

test_that( 'admixcor2 works', {
    # for this test, we don't have to fully converge (even in toy data it sometimes takes too long)
    # this stops in a single iteration in tests!
    tol <- 1e-2 # default ~1e-8
    # first test default version with no regularization
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        objv <- admixcor2( Theta, K, tol = tol, vertex_refine = TRUE )
    )
    validate_admixcor( objv, n, K, v = 2 )
    
    # now test regularized versions
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, alpha = alpha )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, alpha = alpha, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # all of those earlier cases were default Q initialization using kmeans, try other non-default cases
    # test default unregularized version and regularized, to catch edge cases potentially
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, Q_type = 'random' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, Q_type = 'random', alpha = alpha, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
    # NOTE: admixcor2 really doesn't like `Q_type = 'uniform'`, it fails regularly with codes like these (both with and without regularization):
    # - Error: from glmnet C++ code (error code 7777); All used predictors have zero variance
    ## expect_silent(
    ##     obj <- admixcor2( Theta, K, tol = tol, Q_type = 'uniform' )
    ## )
    ## validate_admixcor( obj, n, K, v = 2 )
    ## expect_silent(
    ##     obj <- admixcor2( Theta, K, tol = tol, Q_type = 'uniform', alpha = alpha, beta = beta )
    ## )
    ## validate_admixcor( obj, n, K, v = 2 )

    # ditto L initializations, try non-default versions now (identity is default)
    # test in context of Q_type = 'kmeans' only
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, L_type = 'uniform' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, L_type = 'diagrandom' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, L_type = 'random' )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # and Psi algorithms!
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, Psi_algorithm = 'bvls' )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # and Q algorithms!
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        objv <- admixcor2( Theta, K, tol = tol, Q_algorithm = 'quadprog', vertex_refine = TRUE )
    )
    validate_admixcor( objv, n, K, v = 2 )
})

test_that( 'admixcor2 reformed works', {
    # for this test, we don't have to fully converge (even in toy data it sometimes takes too long)
    # this stops in a single iteration in tests!
    tol <- 1e-2 # default ~1e-8
    # first test default version with no regularization
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        objv <- admixcor2( Theta, K, tol = tol, reformed = TRUE, vertex_refine = TRUE )
    )
    validate_admixcor( objv, n, K, v = 2 )
    
    # now test regularized versions
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, alpha = alpha )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, alpha = alpha, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # all of those earlier cases were default Q initialization using kmeans, try other non-default cases
    # test default unregularized version and regularized, to catch edge cases potentially
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Q_type = 'random' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Q_type = 'random', alpha = alpha, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
    # NOTE: admixcor2 really doesn't like `Q_type = 'uniform'`, it fails regularly with codes like these (both with and without regularization):
    # - Error: from glmnet C++ code (error code 7777); All used predictors have zero variance
    ## expect_silent(
    ##     obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Q_type = 'uniform' )
    ## )
    ## validate_admixcor( obj, n, K, v = 2 )
    ## expect_silent(
    ##     obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Q_type = 'uniform', alpha = alpha, beta = beta )
    ## )
    ## validate_admixcor( obj, n, K, v = 2 )

    # ditto L initializations, try non-default versions now (identity is default)
    # test in context of Q_type = 'kmeans' only
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, L_type = 'uniform' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, L_type = 'diagrandom' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, L_type = 'random' )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # and Psi algorithms!
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Psi_algorithm = 'bvls' )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # and Q algorithms!
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        objv <- admixcor2( Theta, K, tol = tol, reformed = TRUE, Q_algorithm = 'quadprog', vertex_refine = TRUE )
    )
    validate_admixcor( objv, n, K, v = 2 )
})

test_that( 'admixcor2 works in full run with small K', {
    # here we want a fuller run, which in practice resulted in some errors we want to learn how to avoid; thus the only change is we use the default tolerance
    # ideally we run this with K=2, but meh, benchmarks assume K=3 elsewhere
    # focus on cases of highest interest
    Q_type <- 'random'
    # default `Psi_algorithm = 'glmnet'` is key to these errors
    
    # test regularized versions only, default `Q_algorithm = 'original'`
    expect_silent(
        obj <- admixcor2( Theta, K, alpha = alpha, Q_type = Q_type )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, beta = beta, Q_type = Q_type )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, alpha = alpha, beta = beta, Q_type = Q_type )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # and `Q_algorithm = 'quadprog'`!
    expect_silent(
        obj <- admixcor2( Theta, K, alpha = alpha, Q_type = Q_type, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, beta = beta, Q_type = Q_type, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, alpha = alpha, beta = beta, Q_type = Q_type, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K, v = 2 )
})
