# for convenience, take objective in same format as quadprog
# no need to specify constraints, our vertices always satisfy them!
# NOTE: this solves problem for a single individual
update_Q_vertex_scan <- function( D, d, fast = TRUE ) {
    # infer K from inputs
    K <- nrow( D )

    if ( fast ) {
        # can calculate all objectives quickly for these special cases
        objs <- diag( D ) / 2 - d
        # identify the minimum
        i <- which.min( objs )
        # construct the solution
        # this makes a q vector that is zero everywhere except 1 at i
        q_best <- rep.int( 0, K )
        q_best[i] <- 1
        obj_best <- objs[i]
    } else {
        
        # initialize best solution so far (guaranteed to be beat immediately)
        q_best <- NA
        obj_best <- Inf
        
        # navigate vertices
        for ( i in 1 : K ) {
            # this makes a q vector that is zero everywhere except 1 at i
            q <- rep.int( 0, K )
            q[i] <- 1
            # see its objective value
            obj <- obj_quadprog( D, d, q )
            # compare to previous best
            if ( obj < obj_best ) {
                obj_best <- obj
                q_best <- q
            }
        }
    }
    # return best vertex and its objective value
    return( list( q = q_best, obj = obj_best ) )
}

obj_quadprog <- function( D, d, q )
    drop( q %*% D %*% q / 2 - d %*% q )
