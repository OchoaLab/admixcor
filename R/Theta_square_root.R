
Theta_square_root<-function(Theta,K)
{
	ev<-eigen(Theta)
	evec<-ev$vectors[,1:K]
	eval<-diag(sqrt(ev$values[1:K]))
	ThetaSR<-evec%*%eval
	return(ThetaSR)
}
