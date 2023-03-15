admixcor <- function( Theta, K, gamma = 0.01, stop = 1e-15, verbose = FALSE ) {
    n<-nrow(Theta)
    normz<-norm(Theta,"F")

    ThetaSR<-Theta_square_root(Theta,K)
    Vars<-Initialize(Theta,K,n)

    I<-diag(1,K,K)
    
    ndQ<-100
    ndPsiSR<-100
    ndR<-100

    Q0<-Vars[[1]]
    PsiSR0<-Vars[[2]]
    R0<-Vars[[3]]
    f0<-Vars[[4]]
    minf<-100

    nstep<-0
    while(ndQ>stop && ndPsiSR>stop && ndR>stop)
    {
	Qinv<-MASS::ginv(Q0)
	R1<-MASS::ginv(PsiSR0)%*%Qinv%*%ThetaSR
	svd<-svd(R1)
        u<-svd$u
        v<-svd$v
        R1<-u%*%t(v)

	QT<-t(Q0)
	PsiSR1<-MASS::ginv(QT%*%Q0+gamma*I)%*%QT%*%ThetaSR%*%t(R1)
	#PsiSR1<-MASS::ginv(Q1)%*%ThetaSR%*%t(R1)
	PsiSR1[lower.tri(PsiSR1)] <- 0
	PsiSR1[PsiSR1<0]<-0
	PsiSR1[PsiSR1>1]<-1
	
	Q1<-ThetaSR%*%t(R1)%*%MASS::ginv(PsiSR1)
	for(j in 1:n){
            #Q1[j,]=ifelse(Q1[j,]<0.0,Q1[j,]-min(Q1[j,])+0.0001,Q1[j,])
            Q1[j,]<-projsplx(Q1[j,])
	}
	#Q1<-Q1/rowSums(Q1)

	ndQ<-norm((Q0-Q1),"F")/sqrt(n*K)
	ndPsiSR<-norm((PsiSR0-PsiSR1),"F")/K
	ndR<-norm((R0-R1),"F")/K
	f1<-Objective(Theta, Q1, PsiSR1, n, K)

	Q0<-Q1
	PsiSR0<-PsiSR1
	R0<-R1
	f0<-f1

	if(nstep%%1000==0) 
	{
            if (verbose)
                message(ndQ,' ',ndPsiSR,' ',ndR,' ',f1)
            #print(Q0)
            #print(PsiSR0)
	}
	nstep<-nstep+1
	if(minf>f0)
	{
            minf<-f0
            Qfinal<-Q0
            PsiSRfinal<-PsiSR0
	}
	if(nstep>100000) break

    }

    Psifinal<-PsiSRfinal%*%t(PsiSRfinal)
    return(
        list(
            Q = Qfinal,
            Psi = Psifinal
        )
    )
}
