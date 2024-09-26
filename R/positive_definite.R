#Check if Psi matrix is positive definite
#If not modify the Psi matrix to be positive definite

# reg = TRUE applies a common regularization across eigenvalues, which could be overkill depending on how big the difference is between the minimum eigenvalue and the tolerance, and ends up affecting all other higher eigenvalues
# reg = FALSE only edits the offending eigenvalues, leaves the rest intact
positive_definite <- function( Psi, tol_psi = sqrt( .Machine$double.eps ), reg = FALSE ){
    if ( reg ) {
        # for now only values are used, in that case this is faster
        evd <- eigen( Psi, only.values = TRUE )
        min_ev <- min( evd$values )
        
        if ( min_ev <= tol_psi ) {
            K <- nrow( Psi )
            I <- diag( 1, K, K )
            k <- tol_psi - min_ev
            Psi <- Psi + k*I
        }
    } else {
        # use full eigendecomposition here
        evd <- eigen( Psi )
        eva <- evd$values
        # identify eigenvalues that are too small
        indexes <- eva < tol_psi
        # edit only if needed
        if ( any( indexes ) ) {
            # replace those only with our minimum
            eva[ indexes ] <- tol_psi
            # reform Psi with these better eigenvalues
            Psi <- tcrossprod( evd$vectors %*% diag( eva ), evd$vectors )
        }
    }

    return( Psi )
}

