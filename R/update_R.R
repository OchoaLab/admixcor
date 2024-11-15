# updates R using the linear regression formulas, then projects to rotation space
update_R <- function( ThetaSR, Q, L ) {
    # update R (rotation)
    Qinv <- MASS::ginv(Q) # general pseudoinverse, was fastest
    R <- MASS::ginv(L) %*% Qinv %*% ThetaSR # general pseudoinverse, was fastest
    
    # project updated R to rotation space
    # simply replaces eigenvalues with ones
    svd <- svd( R )
    R <- tcrossprod( svd$u, svd$v )
    
    return( R )
}
