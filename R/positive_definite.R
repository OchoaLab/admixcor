#Check if Psi matrix is positive definite
#If not modify the Psi matrix to be positive definite

positive_definite <- function(Psi, K, tol_psi){
	ev <- eigen(Psi)
	
	if(any(ev$values <= 0)) {
		I <- diag( 1, K, K )
		k <- tol_psi - min( ev$values )  
		Psi <- Psi + k*I
	}

	return(Psi)
}
