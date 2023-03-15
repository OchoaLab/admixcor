Initialize<-function(Theta,ThetaSR,k,n)
{
	#Q<-matrix(runif(n*k),ncol=k,nrow=n)
	#Q<-Q/rowSums(Q)

	Q<-matrix(0,nrow=n,ncol=k)
        fit<-kmeans(Theta,k,nstart=100,iter.max=100)
        for(x in 1:n){
	        for(y in 1:k){
		        if(fit$cluster[x]==y) {
			        Q[x,y]<-1
			}
		}
	}

	#for(i in 1:n)
        #{
	#	Q[i,1]<-0.5
	#	Q[i,2]<-0.5
	#}

	#PsiSR<-diag(runif(k),k,k)
	#PsiSR<-matrix(runif(k*k),k,k)
	PsiSR<-matrix(1,k,k)
	#PsiSR<-diag(1,k,k)
	PsiSR[lower.tri(PsiSR)] <- 0
	
	R<-diag(1,k,k)

	f<-Objective(Theta,Q,PsiSR,n,K)

	result<-list(Q=Q,PsiSR=PsiSR,R=R,f=f)

	return(result)
}

