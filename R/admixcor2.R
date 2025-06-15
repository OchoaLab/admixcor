#' @export
admixcor2 <- function(
                     Theta,
                     K,
		     alpha = 0,
                     gamma = 0,
                     L_type = c('diagrandom', 'random'),
                     tol = sqrt( .Machine$double.eps ), # 1e-15
		     tol_psi = sqrt( .Machine$double.eps ), # 1e-15
                     nstep_max = 100000,
                     report_freq = 1000,
		     reformed = FALSE,
		     ties_none = FALSE,
		     tol_stretch = -0.01
                     ) {
    # process options
    L_type <- match.arg( L_type )
    
    # number of individuals is dimensions of Theta
    n <- nrow( Theta )

    # decompose Theta for "linear" objective
    ThetaSR <- theta_square_root( Theta, K )
    # normalize gamma to be on a more similar scale to delta, etc
    # takes into account that dim(L) = (K, K) is much smaller than dim(Q) = dim(ThetaSR) = (n, K)
    gamma <- gamma * n / K
    alpha <- alpha * n / K

    # initialize other variables
    Vars <- initialize( ThetaSR, K, n, L_type )
    Q0 <- Vars$Q
    L0 <- Vars$L
    Psi0 <- tcrossprod( L0 )
    R0 <- Vars$R
    
    if(reformed)
	Theta <- ThetaSR %*% t(ThetaSR)

    # calculate objective of this initial (very bad) solution
    f0 <- objective( ThetaSR, Q0, L0, R0, gamma )
    #g0 <- objective2( Theta, Q0, Psi0, alpha )

    # constant used in regularized expressions
    I <- diag( 1, K, K )
    
    # initialize gradient norms to ensure first iteration occurs
    ndQ <- Inf
    f_best <- Inf
    g_best <- Inf

    # initialize counter
    nstep <- 0

    # initialize select values to save for report
    # values of starting point
    nsteps <- nstep
    ndQs <- NA
    objs <- f0
    dobjs <- NA
    # zeroeth iteration tells me if this was initialized as singular!
    L_singulars <- inherits(try(solve(L0), silent = TRUE), "try-error")

    while ( ndQ > tol ) {
        # increment counter, break if needed
        nstep <- nstep + 1
	if ( nstep > nstep_max )
            break
        
        # apply the updates, one at the time
        Psi1 <- update_Psi( Theta, Q0, alpha )
        Psi1 <- positive_definite( Psi1, tol_psi = tol_psi )
        L1 <- t(chol( Psi1 ) )
        R1 <- update_R( ThetaSR, Q0, L1 )
        obj <- update_Q( ThetaSR, L1, R1, I )
        Q1 <- obj$Q
        L_singular1 <- obj$L_singular
        # if we encounter a singular case we scrap the current solution and draw a completely new one, essentially start from scratch again (but without restarting iteration count)
        # while we could re-draw only Q, this risks a bad L still influencing the R we pick and therefore the next Q, so instead we re-draw everything!
        if ( L_singular1 ) {
            Vars <- initialize( ThetaSR, K, n, L_type )
            Q1 <- Vars$Q
            L1 <- Vars$L
            R1 <- Vars$R
            # NOTE: in the log we will know if a restart occurred because the iteration will be marked with `L_singular = TRUE`
        } else {
            # (skip stretching if we reset)
            # apply stretching
            Q1 <- stretch_Q( Q1, ties_none = ties_none, tol = tol_stretch )$Q
            # adjust stretching due to tolerance for negative values
            Q1[Q1 < 0] <- 0
            Q1 <- Q1 / rowSums( Q1 )
        }
        
        # calculate step sizes, to assess convergence
	ndQ <- norm( Q0 - Q1, "F" )^2
        
        # calculate new objective
	f1 <- objective( ThetaSR, Q1, L1, R1, gamma )
        # and its delta too (this matches norm formulations above)
        # use first value only (total objective, including penalty terms!)
        df <- abs( f1[1L] - f0[1L] )

        # update variables for next iteration
	Q0 <- Q1
	Psi0 <- Psi1
	R0 <- R1
	f0 <- f1

        # update report infrequently
	if ( nstep %% report_freq == 0) {
            # increment report
            nsteps <- c( nsteps, nstep )
            ndQs <- c( ndQs, ndQ )
            objs <- rbind( objs, f0 ) # this is a matrix!
            dobjs <- c( dobjs, df )
            L_singulars <- c( L_singulars, L_singular1 )
	}

        # decide if this is better than previous best, based on total objective (including penalties)
	if ( f_best[1L] > f0[1L] ) {
            f_best <- f0
            Q_best <- Q0
            Psi_best <- Psi0
            R_best <- R0
            # nice extra info, mostly for troubleshooting
            nstep_best <- nstep
            ndQ_best <- ndQ
            df_best <- df
            L_singular_best <- L_singular1
	}
    }

    # report final iteration info
    nsteps <- c( nsteps, nstep )
    ndQs <- c( ndQs, ndQ )
    objs <- rbind( objs, f0 ) # this is a matrix!
    dobjs <- c( dobjs, df )
    L_singulars <- c( L_singulars, L_singular1 )
    # compare data for best
    nsteps <- c( nsteps, nstep_best )
    ndQs <- c( ndQs, ndQ_best )
    objs <- rbind( objs, f_best ) # this is a matrix!
    dobjs <- c( dobjs, df_best )
    L_singulars <- c( L_singulars, L_singular_best )
    # finalize report
    report <- tibble::tibble(
                          nstep = nsteps,
                          dQ = ndQs,
                          # treat matrix cases accordingly, extracting columns as needed
                          obj = objs[,1L],
                          dobj = dobjs,
                          O = objs[,2L],
                          pPsi = objs[,3L],
                          L_singular = L_singulars
                      )
    
    
    return(
        list(
            Q = Q_best,
            Psi = Psi_best,
            f = f_best[1L], # only first value makes most sense in this context
            report = report,
            # for ease of troubleshooting, all internal values needed to calculate on the outside the internal objective function
            ThetaSR = ThetaSR,
            R = R_best
        )
    )
}
