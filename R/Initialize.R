# initializes Q using kmeans on Theta
# L initialized to fixed matrix of 1/k's upper triangular
# R is initialized to I
Initialize<-function( Theta, ThetaSR, K, n, gamma, delta ) {
    # initialize Q randomly
    #Q<-matrix(runif(n*k),ncol=k,nrow=n)
    #Q<-Q/rowSums(Q)

    # initialize Q with k-means!
    Q<-matrix(0,nrow=n,ncol=K)
    fit<-stats::kmeans(Theta,K,nstart=100,iter.max=100)
    for(x in 1:n){
        for(y in 1:K){
            if(fit$cluster[x]==y) {
                Q[x,y]<-1
            }
        }
    }

    # initialize L (various versions)
    #L<-diag(runif(K),K,K)
    #L<-matrix(runif(K*K),K,K)
    L <- matrix( 1/K, K, K ) # make sure final Psi values don't exceed 1
    #L<-diag(1,K,K)
    L[lower.tri(L)] <- 0

    # initialize R to identity matrix
    R <- diag( 1, K, K )

    # calculate objective of this initial (very bad) solution
    f <- Objective( ThetaSR, Q, L, R, gamma, delta )

    # return everything!
    return(
        list(
            Q = Q,
            L = L,
            R = R,
            f = f
        )
    )
}

