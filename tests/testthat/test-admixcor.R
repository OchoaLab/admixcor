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
# some default values for the regularization parameters
gamma <- 0.01
delta <- 0.01
# constant used in regularized expressions
I <- diag( 1, K, K )

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
validate_L <- function( L, K ) {
    expect_true( is.matrix( L ) )
    expect_equal( nrow( L ), K )
    expect_equal( ncol( L ), K )
    expect_true( min( L ) >= 0 )
    expect_true( max( L ) <= 1 )
    # demand that lower triangle is zero! (excludes diagonal)
    expect_true( all( L[ lower.tri( L ) ] == 0 ) )
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

# uniform testing for all cases
validate_initialize <- function( obj, n, K ) {
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('Q', 'L', 'R') )
    validate_Q( obj$Q, n, K )
    validate_R( obj$R, K, I )
    validate_L( obj$L, K )
}

# uniform testing for all cases
validate_admixcor <- function( obj, n, K ) {
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('Q', 'Psi', 'f', 'report', 'ThetaSR', 'L', 'R') )
    validate_Q( obj$Q, n, K )
    validate_L( obj$L, K )
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
L2[ lower.tri( L2 ) ] <- 0


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
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use exact ThetaSR from real data, but random Q for test
    L2 <- update_L( ThetaSR, Q2, R, gamma, I )
    # test basic expectations
    validate_L( L2, K )

    # now try exact solution (true ThetaSR, Q and R), here it is recoverable in theory if we don't penalize!
    L2 <- update_L( ThetaSR, Q, R, 0, I )
    expect_equal( L2, L )

    # test new version based on nnls, repeating both earlier tests
    # here it's fine to exclude gamma and I
    L2 <- update_L( ThetaSR, Q2, R, algorithm = 'nnls' )
    validate_L( L2, K )
    L2 <- update_L( ThetaSR, Q, R, algorithm = 'nnls' )
    expect_equal( L2, L )

    # ditto bvls
    L2 <- update_L( ThetaSR, Q2, R, algorithm = 'bvls' )
    validate_L( L2, K )
    L2 <- update_L( ThetaSR, Q, R, algorithm = 'bvls' )
    expect_equal( L2, L )
    
    # ditto glmnet
    L2 <- update_L( ThetaSR, Q2, R, gamma, I, algorithm = 'glmnet' )
    validate_L( L2, K )
    # can only recover true solution if we don't penalize!
    L2 <- update_L( ThetaSR, Q, R, 0, I, algorithm = 'glmnet' )
    expect_equal( L2, L )
})

test_that( 'update_Q works', {
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use exact ThetaSR from real data, but random L for test
    Q2 <- update_Q( ThetaSR, L2, R, delta, I )
    # test basic expectations
    validate_Q( Q2, n, K )

    # now try exact solution (true ThetaSR, L and R), here it is recoverable in theory if we don't penalize!
    Q2 <- update_Q( ThetaSR, L, R, 0, I )
    expect_equal( Q2, Q )

    # test new version based on nnls, repeating both earlier tests
    # here it's fine to exclude delta and I
    Q2 <- update_Q( ThetaSR, L2, R, algorithm = 'nnls' )
    validate_Q( Q2, n, K )
    Q2 <- update_Q( ThetaSR, L, R, algorithm = 'nnls' )
    expect_equal( Q2, Q )

    # ditto bvls
    Q2 <- update_Q( ThetaSR, L2, R, algorithm = 'bvls' )
    validate_Q( Q2, n, K )
    Q2 <- update_Q( ThetaSR, L, R, algorithm = 'bvls' )
    expect_equal( Q2, Q )

    # ditto glmnet
    Q2 <- update_Q( ThetaSR, L2, R, delta, I, algorithm = 'glmnet' )
    validate_Q( Q2, n, K )
    # can only recover true solution if we don't penalize!
    Q2 <- update_Q( ThetaSR, L, R, 0, I, algorithm = 'glmnet' )
    expect_equal( Q2, Q )

    # ditto quadprog
    Q2 <- update_Q( ThetaSR, L2, R, delta, I, algorithm = 'quadprog' )
    validate_Q( Q2, n, K )
    # can only recover true solution if we don't penalize!
    Q2 <- update_Q( ThetaSR, L, R, 0, I, algorithm = 'quadprog' )
    expect_equal( Q2, Q )

    # ditto compact version of quadprog
    Q2 <- update_Q( ThetaSR, L2, R, delta, I, algorithm = 'quadprog-compact' )
    validate_Q( Q2, n, K )
    # can only recover true solution if we don't penalize!
    Q2 <- update_Q( ThetaSR, L, R, 0, I, algorithm = 'quadprog-compact' )
    expect_equal( Q2, Q )
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
    
    # now test no regularization version, where some things may be singular if we're not careful
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = 0 )
    )
    validate_admixcor( obj, n, K )

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
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_type = 'uniform', gamma = 0, delta = 0 )
    )
    validate_admixcor( obj, n, K )

    # ditto L initializations, try non-default versions now (identity is default)
    # test in context of Q_type = 'kmeans' only
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'uniform' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'diagrandom' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'random' )
    )
    validate_admixcor( obj, n, K )

    # and L algorithms!
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_algorithm = 'nnls' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_algorithm = 'bvls' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_algorithm = 'glmnet' )
    )
    validate_admixcor( obj, n, K )

    # and Q algorithms!
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'nnls' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'bvls' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'glmnet' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'quadprog' )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, Q_algorithm = 'quadprog-compact' )
    )
    validate_admixcor( obj, n, K )
})
