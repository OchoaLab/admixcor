

Objective<-function(ThetaSR, Q, PsiSR, R)
{
	return(norm(ThetaSR-(Q%*%PsiSR%*%R)))
}
