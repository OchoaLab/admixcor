

Objective<-function(Theta, Q, PsiSR, n, K)
{
	Psi<-PsiSR%*%t(PsiSR)
	return(norm(Theta-(Q%*%Psi%*%t(Q)))/sqrt(n*K))
}
