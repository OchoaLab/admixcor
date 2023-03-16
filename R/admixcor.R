admixcor <- function(
                     Theta,
                     K,
                     gamma = 0.01,
                     delta = 0.01,
                     stop = 1e-15,
                     nstep_max = 100000,
                     report_freq = 1000
                     ) {
    n<-nrow(Theta)

    # decompose Theta for "linear" objective
    ThetaSR<-Theta_square_root(Theta,K)
    # normalization factor that puts penalty terms in equal footing
    normz <- norm( ThetaSR, "F" )^2
    # apply to parameters
    gamma <- gamma * normz
    delta <- delta * normz
    #message( 'normz: ', normz )
    
    # initialize other variables
    Vars <- Initialize( Theta, K, n )
    Q0 <- Vars$Q
    L0 <- Vars$L
    R0 <- Vars$R

    # calculate objective of this initial (very bad) solution
    f0 <- Objective( ThetaSR, Q0, L0, R0, gamma, delta )

    # constant used in regularized expressions
    I<-diag(1,K,K)

    # initialize gradient norms to ensure first iteration occurs
    ndQ <- Inf
    ndL <- Inf
    ndR <- Inf
    f_best <- Inf

    # initialize counter
    nstep<-0

    # initialize select values to save for report
    # values of starting point
    nsteps <- -1
    ndQs <- NA
    ndLs <- NA
    ndRs <- NA
    objs <- f0
    dobjs <- NA

    while (ndQ>stop && ndL>stop && ndR>stop) {
        # apply the updates, one at the time
        R1 <- update_R( ThetaSR, Q0, L0 )
        L1 <- update_L( ThetaSR, Q0, R1, gamma, I )
        Q1 <- update_Q( ThetaSR, L1, R1, delta, I )
        
        # calculate step sizes, to assess convergence
	ndQ<-norm((Q0-Q1),"F")/sqrt(n*K)
	ndL<-norm((L0-L1),"F")/K
	ndR<-norm((R0-R1),"F")/K
        
        # calculate new objective
	f1 <- Objective( ThetaSR, Q1, L1, R1, gamma, delta )
        # and its delta too (this matches norm formulations above)
        # use first value only (total objective, including penalty terms!)
        df <- abs( f1[1L] - f0[1L] )

        # update variables for next iteration
	Q0<-Q1
	L0<-L1
	R0<-R1
	f0<-f1

        # update report infrequently
	if ( nstep %% report_freq == 0) {
            # increment report
            nsteps <- c( nsteps, nstep )
            ndQs <- c( ndQs, ndQ )
            ndLs <- c( ndLs, ndL )
            ndRs <- c( ndRs, ndR )
            objs <- rbind( objs, f0 ) # this is a matrix!
            dobjs <- c( dobjs, df )
            #print(Q0)
            #print(L0)
	}

        # decide if this is better than previous best, based on total objective (including penalties)
	if ( f_best[1L] > f0[1L] ) {
            f_best<-f0
            Q_best<-Q0
            L_best<-L0
            # nice extra info, mostly for troubleshooting
            nstep_best <- nstep
            ndQ_best <- ndQ
            ndL_best <- ndL
            ndR_best <- ndR
            df_best <- df
	}
	nstep<-nstep+1
	if (nstep>nstep_max) break
    }

    # report final iteration info
    # nstep-1 because we incremented without a change in last case (and it starts from zero because we don't want to change things)
    # increment report
    nsteps <- c( nsteps, nstep-1 )
    ndQs <- c( ndQs, ndQ )
    ndLs <- c( ndLs, ndL )
    ndRs <- c( ndRs, ndR )
    objs <- rbind( objs, f0 ) # this is a matrix!
    dobjs <- c( dobjs, df )
    # compare data for best
    nsteps <- c( nsteps, nstep_best )
    ndQs <- c( ndQs, ndQ_best )
    ndLs <- c( ndLs, ndL_best )
    ndRs <- c( ndRs, ndR_best )
    objs <- rbind( objs, f_best ) # this is a matrix!
    dobjs <- c( dobjs, df_best )
    # finalize report
    report <- tibble::tibble(
                          nstep = nsteps,
                          dQ = ndQs,
                          dL = ndLs,
                          dR = ndRs,
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
