admixcor <- function(
                     Theta,
                     K,
                     gamma = 0.01,
                     delta = 0.01,
                     Q_type = c('kmeans', 'random', 'uniform'),
                     L_type = c('identity', 'uniform', 'diagrandom', 'random'),
                     tol = sqrt( .Machine$double.eps ), # 1e-15
                     nstep_max = 100000,
                     report_freq = 1000
                     ) {
    # process options
    Q_type <- match.arg( Q_type )
    L_type <- match.arg( L_type )

    # number of individuals is dimensions of Theta
    n <- nrow( Theta )

    # decompose Theta for "linear" objective
    ThetaSR <- Theta_square_root( Theta, K )
    # normalize gamma to be on a more similar scale to delta, etc
    # takes into account that dim(L) = (K, K) is much smaller than dim(Q) = dim(ThetaSR) = (n, K)
    gamma <- gamma * n / K
    
    # initialize other variables
    Vars <- Initialize( Theta, K, n, Q_type, L_type )
    Q0 <- Vars$Q
    L0 <- Vars$L
    R0 <- Vars$R

    # calculate objective of this initial (very bad) solution
    f0 <- Objective( ThetaSR, Q0, L0, R0, gamma, delta )

    # constant used in regularized expressions
    I <- diag( 1, K, K )
    
    # initialize gradient norms to ensure first iteration occurs
    ndQ <- Inf
    f_best <- Inf

    # initialize counter
    nstep <- 0

    # initialize select values to save for report
    # values of starting point
    nsteps <- nstep
    ndQs <- NA
    objs <- f0
    dobjs <- NA

    while ( ndQ > tol ) {
        # increment counter, break if needed
        nstep <- nstep + 1
	if ( nstep > nstep_max )
            break
        
        # apply the updates, one at the time
        R1 <- update_R( ThetaSR, Q0, L0 )
        L1 <- update_L( ThetaSR, Q0, R1, gamma, I )
        Q1 <- update_Q( ThetaSR, L1, R1, delta, I )
        
        # calculate step sizes, to assess convergence
	ndQ <- norm( Q0 - Q1, "F" )^2
        
        # calculate new objective
	f1 <- Objective( ThetaSR, Q1, L1, R1, gamma, delta )
        # and its delta too (this matches norm formulations above)
        # use first value only (total objective, including penalty terms!)
        df <- abs( f1[1L] - f0[1L] )

        # update variables for next iteration
	Q0 <- Q1
	L0 <- L1
	R0 <- R1
	f0 <- f1

        # update report infrequently
	if ( nstep %% report_freq == 0) {
            # increment report
            nsteps <- c( nsteps, nstep )
            ndQs <- c( ndQs, ndQ )
            objs <- rbind( objs, f0 ) # this is a matrix!
            dobjs <- c( dobjs, df )
	}

        # decide if this is better than previous best, based on total objective (including penalties)
	if ( f_best[1L] > f0[1L] ) {
            f_best <- f0
            Q_best <- Q0
            L_best <- L0
            # nice extra info, mostly for troubleshooting
            nstep_best <- nstep
            ndQ_best <- ndQ
            df_best <- df
	}
    }

    # report final iteration info
    nsteps <- c( nsteps, nstep )
    ndQs <- c( ndQs, ndQ )
    objs <- rbind( objs, f0 ) # this is a matrix!
    dobjs <- c( dobjs, df )
    # compare data for best
    nsteps <- c( nsteps, nstep_best )
    ndQs <- c( ndQs, ndQ_best )
    objs <- rbind( objs, f_best ) # this is a matrix!
    dobjs <- c( dobjs, df_best )
    # finalize report
    report <- tibble::tibble(
                          nstep = nsteps,
                          dQ = ndQs,
                          # treat matrix cases accordingly, extracting columns as needed
                          obj = objs[,1L],
                          dobj = dobjs,
                          O = objs[,2L],
                          pL = objs[,3L],
                          pQ = objs[,4L]
                      )
    
    # compose final Psi
    Psi_best <- tcrossprod( L_best )
    
    return(
        list(
            Q = Q_best,
            Psi = Psi_best,
            f = f_best[1L], # only first value makes most sense in this context
            report = report
        )
    )
}
