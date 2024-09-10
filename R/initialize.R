# initializes Q using kmeans on Theta (will also work and is faster if ThetaSR is provided!)
# L initialized to I
# R is initialized to I

# NOTES:
# - Q_type:
#   - 'kmeans' (default) continues to be best in latest tests
#   - 'random' is not bad but usually a little worse than 'kmeans'
#   - 'uniform' is super bad, always converged in one iteration (current operation order and other options, perhaps this changes), so it appears to be fixed point.  Will leave in order to run more tests later
initialize <- function(
                       Theta,
                       K,
                       n,
                       Q_type = c('kmeans', 'random', 'uniform'),
                       L_type = c('identity', 'uniform', 'diagrandom', 'random'),
                       fix_L = FALSE
                       ) {
    # process options
    Q_type <- match.arg( Q_type )
    L_type <- match.arg( L_type )

    # initialize Q
    if ( Q_type == 'kmeans' ) {
        # initialize Q with k-means!
        clusters <- stats::kmeans( Theta, K, nstart = 100, iter.max = 100 )
        clusters <- clusters$cluster
        # translate cluster vector to Q matrix
        Q <- matrix( 0, nrow = n, ncol = K )
        for ( x in 1L : n ) {
            Q[ x, clusters[x] ] <- 1
        }
    } else if ( Q_type == 'random' ) {
        # initialize Q randomly
        Q <- matrix( stats::runif( n * K ), ncol = K, nrow = n )
        Q <- Q / rowSums( Q )
    } else if ( Q_type == 'uniform' ) {
        # start in the middle of the simplex for every individual
        # NOTE this is a singular case!
        Q <- matrix( 1 / K, ncol = K, nrow = n )
    } else
        stop( '`Q_type` not implemented: ', Q_type )

    # initialize L
    if ( L_type == 'uniform' ) {
        # make sure final Psi values don't exceed 1
        # all values 1/sqrt(K), ensure Psi is between 0 and 1
        L <- matrix( 1 / sqrt( K ), K, K )
        # ensure this is like cholesky
        if ( fix_L ) {
            L[ upper.tri( L ) ] <- 0
        } else {
            L[ lower.tri( L ) ] <- 0
        }
    } else if ( L_type == 'identity' ) {
        L <- diag( 1, K, K )
    } else if ( L_type == 'diagrandom' ) {
        L <- diag( stats::runif( K ), K, K )
    } else if ( L_type == 'random' ) {
        # random values between 0 and 1/sqrt(K), ensure Psi is between 0 and 1
        L <- matrix( stats::runif( K^2 ), K, K ) / sqrt(K)
        # ensure this is like cholesky
        if ( fix_L ) {
            L[ upper.tri( L ) ] <- 0
        } else {
            L[ lower.tri( L ) ] <- 0
        }
    } else
        stop( '`L_type` not implemented: ', L_type )
    
    # initialize R to identity matrix
    R <- diag( 1, K, K )

    # return everything!
    return(
        list(
            Q = Q,
            L = L,
            R = R
        )
    )
}

