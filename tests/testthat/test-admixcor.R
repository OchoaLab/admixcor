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

test_that( 'Theta_square_root works', {
    # does a dimensionality reduction too, so it's not an exact square root unless it was low rank and it is correct
    expect_silent(
        Theta_sqrt <- Theta_square_root( Theta, K )
    )
    expect_true( is.matrix( Theta_sqrt ) )
    expect_equal( nrow( Theta_sqrt ), n )
    expect_equal( ncol( Theta_sqrt ), K )
    # reconstruct theta, here it will be exact!
    Theta_redone <- tcrossprod( Theta_sqrt )
    expect_equal( Theta_redone, Theta )
})

test_that( 'Objective works', {
    # in this test match is exact, so objective should be zero
    PsiSR <- sqrt( Psi )
    expect_silent( 
        f <- Objective( Theta, Q, PsiSR, n, K )
    )
    expect_equal( f, 0 )

    # now try an inexact match
    PsiSR <- diag( 1, K )
    # passed Q2 instead of Q here too, so it's totally random
    expect_silent( 
        f <- Objective( Theta, Q2, PsiSR, n, K )
    )
    expect_true( is.numeric( f ) )
    expect_equal( length( f ), 1L )
    expect_true( f >= 0 )
    expect_true( f <= 1 ) # holds because of the way it's normalized
})

test_that( 'Initialize works', {
    expect_silent(
        obj <- Initialize( Theta, K, n )
    )
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('Q', 'PsiSR', 'R', 'f') )
    expect_true( is.matrix( obj$Q ) )
    expect_equal( nrow( obj$Q ), n )
    expect_equal( ncol( obj$Q ), K )
    expect_true( min( obj$Q ) >= 0 )
    expect_equal( rowSums( obj$Q ), rep.int( 1, n ) )
    # current initialization results in fixed matrices for PsiSR and R, but better to test for any broadly appropriate values
    expect_true( is.matrix( obj$PsiSR ) )
    expect_equal( nrow( obj$PsiSR ), K )
    expect_equal( ncol( obj$PsiSR ), K )
    expect_true( is.matrix( obj$R ) )
    expect_equal( nrow( obj$R ), K )
    expect_equal( ncol( obj$R ), K )
    expect_true( is.numeric( obj$f ) )
    expect_equal( length( obj$f ), 1L )
    expect_true( obj$f >= 0 )
    expect_true( obj$f <= 1 ) # holds because of the way it's normalized
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
    expect_true( is.matrix( Q3 ) )
    expect_equal( nrow( Q3 ), n )
    expect_equal( ncol( Q3 ), K )
    expect_true( min( Q3 ) >= 0 )
    expect_equal( rowSums( Q3 ), rep.int( 1, n ) )
})

test_that( 'admixcor works', {
    # for this test, we don't have to fully converge (even in toy data it sometimes takes too long)
    # this stops in a single iteration in tests!
    stop <- 1e-2 # default 1e-15
    expect_silent(
        obj <- admixcor( Theta, K, stop = stop )
    )
    ## obj <- admixcor( Theta, K, stop = stop, verbose = TRUE ) # for testing
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('Q', 'Psi') )
    expect_true( is.matrix( obj$Q ) )
    expect_equal( nrow( obj$Q ), n )
    expect_equal( ncol( obj$Q ), K )
    expect_true( min( obj$Q ) >= 0 )
    expect_equal( rowSums( obj$Q ), rep.int( 1, n ) )
    expect_true( is.matrix( obj$Psi ) )
    expect_equal( nrow( obj$Psi ), K )
    expect_equal( ncol( obj$Psi ), K )
    expect_true( min( obj$Psi ) >= 0 )
    # expect_true( max( obj$Psi ) <= 1 ) # current algorithm doesn't enforce this
})
