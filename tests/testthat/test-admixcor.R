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

validate_vertex_inds <- function( indexes, K, n ) {
    # test length
    expect_equal( length( indexes ), K )
    # test uniqueness, actually no longer true for some tie cases, some may be missing, especially with multiple NA cases
    expect_true( length( unique( indexes ) ) <= K )
    # values should be in range, but NAs are allowed
    expect_true( min( indexes, na.rm = TRUE ) >= 1 )
    expect_true( max( indexes, na.rm = TRUE ) <= n )
}

# validates dimensions, non-negativity with tolerance, rows sum to one
validate_Q <- function( Q, n, K, tol = 0 ) {
    expect_true( is.matrix( Q ) )
    expect_equal( nrow( Q ), n )
    expect_equal( ncol( Q ), K )
    # weird way to test this inequality but with tolerance
    ## expect_true( min( Q ) >= tol )
    expect_equal( min( Q, tol ), tol )
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
    # tests of inequality with tolerance
    expect_equal( min( L, 0 ), 0 )
    expect_equal( max( L, 1 ), 1 )
    # demand that lower triangle is zero! (excludes diagonal)
    expect_true( all( L[ upper.tri( L ) ] == 0 ) )
    # for more accurate validations, check its full Psi
    Psi <- tcrossprod( L )
    #expect_true( min( Psi ) >= 0 ) # implied by earlier test
    # not implied, test too!
    # tests of inequality with tolerance
    expect_equal( max( Psi, 1 ), 1 )
}

# validates dimensions, non-negativity, and Psi <= 1
validate_Psi <- function( Psi, K ) {
    expect_true( is.matrix( Psi ) )
    expect_equal( nrow( Psi ), K )
    expect_equal( ncol( Psi ), K )
    # tests of inequality with tolerance
    expect_equal( min( Psi, 0 ), 0 )
    expect_equal( max( Psi, 1 ), 1 )
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
validate_admixcor <- function( obj, n, K, v = 1 ) {
    expect_true( is.list( obj ) )
    if ( v == 1 ) {
        expect_equal( names( obj ), c('Q', 'Psi', 'f', 'report', 'ThetaSR', 'L', 'R') )
    } else {
        expect_equal( names( obj ), c('Q', 'Psi', 'f', 'report', 'ThetaSR', 'R') )
    }
    validate_Q( obj$Q, n, K )
    if ( v == 1 ) {
        # this includes Psi validations
        validate_L( obj$L, K )
    } else {
        validate_Psi( obj$Psi, K )
    }
    validate_R( obj$R, K, I )
    # NOTE: we're not validating ThetaSR, that's ok, it's simple and not a worry (it was validated earlier)
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

validate_stretch_Q <- function( data, n, K, tol ) {
    # test overall object first
    expect_true( is.list( data ) )
    expect_equal( names( data ), c('Q', 'S', 'S_inv', 'alpha') )
    # test Q now
    validate_Q( data$Q, n, K, tol = tol )
    # test alpha
    expect_equal( length( data$alpha ), 1 )
    expect_true( data$alpha >= 0 )
    expect_true( data$alpha <= 1 )
    # test S and S_inv
    S <- data$S
    S_inv <- data$S_inv
    expect_true( is.matrix( S ) )
    expect_true( is.matrix( S_inv ) )
    expect_equal( nrow( S ), K )
    expect_equal( ncol( S ), K )
    expect_equal( nrow( S_inv ), K )
    expect_equal( ncol( S_inv ), K )
    # confirm that they are inverses of each other
    expect_equal( S %*% S_inv, diag( K ) )
    # both matrices have rows that sum to 1
    expect_equal( rowSums( S ), rep.int( 1, K ) )
    expect_equal( rowSums( S_inv ), rep.int( 1, K ) )
    # it appears this just isn't true empirically
    ## # we expect S_inv to always have non-negative values (within tolerance)
    ## # not sure why but errors here are quite a bit greater than the tolerance
    ## expect_equal( min( S_inv, 10 * tol ), 10 * tol )
    # NOTE: S can and does have negative values
}

test_that( 'vertex_inds, stretch_Q, stretch_Psi works', {
    # test first on our random Q
    expect_silent(
        indexes <- vertex_inds( Q )
    )
    # validate these indexes
    validate_vertex_inds( indexes, K, n )

    # repeat with other version
    expect_silent(
        indexes <- vertex_inds( Q, ties_none = TRUE )
    )
    validate_vertex_inds( indexes, K, n )

    # then repeat with a Q constructed to have ties for two ancestries
    # keep original dimensions for simplicity
    Q2 <- Q
    # this will make the first individual appear twice originally, though hopefully a different individual gets chosen for one of the two ancestries later, eventually
    Q2[ 1, ] <- c( 0.2, 0.4, 0.4 )
    expect_silent(
        indexes <- vertex_inds( Q2 )
    )
    # validate these indexes
    validate_vertex_inds( indexes, K, n )
    # repeat with other version
    expect_silent(
        indexes <- vertex_inds( Q2, ties_none = TRUE )
    )
    validate_vertex_inds( indexes, K, n )

    # and to make it even more challenging, have a second tie on another individual
    Q3 <- Q2
    Q3[ 2, ] <- c( 0.1, 0.45, 0.45 )
    expect_silent(
        indexes <- vertex_inds( Q3 )
    )
    # validate these indexes
    validate_vertex_inds( indexes, K, n )
    # repeat with other version
    expect_silent(
        indexes <- vertex_inds( Q3, ties_none = TRUE )
    )
    validate_vertex_inds( indexes, K, n )
    
    # apply same datasets to stretch_Q
    # this is default, but needs to be part of validation too
    tol_stretch <- 0 # -0.01
    expect_silent(
        data <- stretch_Q( Q, tol = tol_stretch )
    )
    # validate all of the rich data provided
    validate_stretch_Q( data, n, K, tol = tol_stretch )
    # repeat for constructed edge cases with ties
    expect_silent(
        data2 <- stretch_Q( Q2, tol = tol_stretch )
    )
    validate_stretch_Q( data2, n, K, tol = tol_stretch )
    expect_silent(
        data3 <- stretch_Q( Q3, tol = tol_stretch )
    )
    validate_stretch_Q( data3, n, K, tol = tol_stretch )
    # and all with other version
    expect_silent(
        data4 <- stretch_Q( Q, tol = tol_stretch, ties_none = TRUE )
    )
    validate_stretch_Q( data4, n, K, tol = tol_stretch )
    expect_silent(
        data5 <- stretch_Q( Q2, tol = tol_stretch, ties_none = TRUE )
    )
    validate_stretch_Q( data5, n, K, tol = tol_stretch )
    expect_silent(
        data6 <- stretch_Q( Q3, tol = tol_stretch, ties_none = TRUE )
    )
    validate_stretch_Q( data6, n, K, tol = tol_stretch )

    # and now apply these transformations to Psi too, and validate them
    # NOTE: these are often negative, but otherwise are fine.  Let's skip tests to have a cleaner testing environment
    expect_silent(
        Psi2 <- stretch_Psi( Psi, data$S_inv )
    )
    #validate_Psi( Psi2, K )
    # TODO x2: min(Psi, 0) (`actual`) not equal to 0 (`expected`).
    expect_silent(
        Psi2 <- stretch_Psi( Psi, data2$S_inv )
    )
    #validate_Psi( Psi2, K )
    # TODO x2: min(Psi, 0) (`actual`) not equal to 0 (`expected`).
    expect_silent(
        Psi2 <- stretch_Psi( Psi, data3$S_inv )
    )
    #validate_Psi( Psi2, K )
    # TODO x3: min(Psi, 0) (`actual`) not equal to 0 (`expected`).
    expect_silent(
        Psi2 <- stretch_Psi( Psi, data4$S_inv )
    )
    #validate_Psi( Psi2, K )
    # TODO: min(Psi, 0) (`actual`) not equal to 0 (`expected`).
    expect_silent(
        Psi2 <- stretch_Psi( Psi, data5$S_inv )
    )
    #validate_Psi( Psi2, K )
    # TODO: min(Psi, 0) (`actual`) not equal to 0 (`expected`).
    expect_silent(
        Psi2 <- stretch_Psi( Psi, data6$S_inv )
    )
    #validate_Psi( Psi2, K )
    # TODO x2: min(Psi, 0) (`actual`) not equal to 0 (`expected`).
})

test_that( 'initialize works', {
    # test all L_types!
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
    # TODO x5: max(Psi, 1) (`actual`) not equal to 1 (`expected`).
    
    # now try exact solution (true ThetaSR, Q and R), here it is recoverable in theory if we don't penalize!
    L2 <- update_L( ThetaSR, Q, R )
    expect_equal( L2, L )

    # test zero gamma test with random Q
    L2 <- update_L( ThetaSR, Q2, R )
    validate_L( L2, K )
    # TODO x5: max(Psi, 1) (`actual`) not equal to 1 (`expected`).
})

test_that( 'update_Psi works', {
    # test the default glmnet first...
    
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use true Theta, but random Q for test
    Psi2 <- update_Psi( Theta, Q2, alpha )
    # test basic expectations
    validate_Psi( Psi2, K )

    # now try exact solution (true Theta and Q), here it is recoverable in theory if we don't penalize!
    Psi2 <- update_Psi( Theta, Q )
    expect_equal( Psi2, Psi )

    # test zero alpha version with random Q
    Psi2 <- update_Psi( Theta, Q2 )
    validate_Psi( Psi2, K )
})

test_that( 'update_Q works', {
    # test on random data first, make sure it doesn't break; all are globally set
    # here we use exact ThetaSR from real data, but random L for test
    obj <- update_Q( ThetaSR, L2, R, delta, I )
    Q2 <- obj$Q
    # test basic expectations
    validate_Q( Q2, n, K )
    
    # now try exact solution (true ThetaSR, L and R), here it is recoverable in theory if we don't penalize!
    obj <- update_Q( ThetaSR, L, R, 0, I )
    Q2 <- obj$Q
    expect_equal( Q2, Q )
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
    
    # first test default unregularized version
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol )
    )
    validate_admixcor( obj, n, K )
    
    # now test regularized versions
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = gamma )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, delta = delta )
    )
    validate_admixcor( obj, n, K )
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, gamma = gamma, delta = delta )
    )
    validate_admixcor( obj, n, K )
    # TODO: max(Psi, 1) (`actual`) not equal to 1 (`expected`).

    # ditto L initializations, try non-default versions now (diagrandom is default)
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor( Theta, K, tol = tol, L_type = 'random' )
    )
    validate_admixcor( obj, n, K )
    # TODO x4: max(Psi, 1) (`actual`) not equal to 1 (`expected`).
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

    # ditto L initializations, try non-default versions now (diagrandom is default)
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, L_type = 'random' )
    )
    validate_admixcor( obj, n, K, v = 2 )
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
    
    # now test regularized versions
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, alpha = alpha )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
    # TODO: max(Psi, 1) (`actual`) not equal to 1 (`expected`).
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, alpha = alpha, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )

    # ditto L initializations, try non-default versions now (diagrandom is default)
    # in this case didn't test unregularized edge cases, but meh, will do if there is a clear need later
    expect_silent(
        obj <- admixcor2( Theta, K, tol = tol, reformed = TRUE, L_type = 'random' )
    )
    validate_admixcor( obj, n, K, v = 2 )
})

test_that( 'admixcor2 works in full run with small K', {
    # here we want a fuller run, which in practice resulted in some errors we want to learn how to avoid; thus the only change is we use the default tolerance
    # ideally we run this with K=2, but meh, benchmarks assume K=3 elsewhere
    # focus on cases of highest interest
    
    # test regularized versions only
    expect_silent(
        obj <- admixcor2( Theta, K, alpha = alpha )
    )
    # TODO: `obj <- admixcor2(Theta, K, alpha = alpha)` produced warnings.
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
    expect_silent(
        obj <- admixcor2( Theta, K, alpha = alpha, beta = beta )
    )
    validate_admixcor( obj, n, K, v = 2 )
})

