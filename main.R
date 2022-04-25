library(MASS)

source('/home/as1052/parameter_fitting/src/linear_model/Theta_square_root.R')
source('/home/as1052/parameter_fitting/src/linear_model/Objective.R')
source('/home/as1052/parameter_fitting/src/linear_model/Initialize.R')
source('/home/as1052/parameter_fitting/src/linear_model/projsplx.R')
source('/home/as1052/parameter_fitting/src/align_Q.R')

args = commandArgs(trailingOnly=TRUE);
if(length(args) !=4 ) {stop("need exact 5 arguements\nRscript AdmixCor kinship.txt(tab delimited file) K(number of populations) ", call.=FALSE); }

coances=args[1]
K=as.numeric(args[2])

Theta<-as.matrix(read.table(coances, header = FALSE, sep = "\t")) #coancestry matrix
Qin<-as.matrix(read.table(args[3], header = FALSE, sep = "\t")) #admixture matrix
Psiin<-as.matrix(read.table(args[4], header = FALSE, sep = "\t")) #Psi matrix
#stop<-as.numeric(args[5])
stop<-1e-15

n<-nrow(Theta)
normz<-norm(Theta,"F")

ptm<-(proc.time())

ThetaSR<-Theta_square_root(Theta,K)
Vars<-Initialize(Theta,ThetaSR,K,n)

I<-diag(1,K,K)
gamma<-0.00

ndQ<-100
ndPsiSR<-100
ndR<-100

Q0<-Vars[[1]]
PsiSR0<-Vars[[2]]
R0<-Vars[[3]]
f0<-Vars[[4]]
minf<-100

ptm <- proc.time()

nstep<-0
while(ndQ>stop && ndPsiSR>stop && ndR>stop)
{
	Qinv<-ginv(Q0)
	R1<-ginv(PsiSR0)%*%Qinv%*%ThetaSR
	svd<-svd(R1)
        u<-svd$u
        v<-svd$v
        R1<-u%*%t(v)

	QT<-t(Q0)
	PsiSR1<-ginv(QT%*%Q0+gamma*I)%*%QT%*%ThetaSR%*%t(R1)
	#PsiSR1<-ginv(Q1)%*%ThetaSR%*%t(R1)
	PsiSR1[lower.tri(PsiSR1)] <- 0
	PsiSR1[PsiSR1<0]<-0
	PsiSR1[PsiSR1>1]<-1
	
	Q1<-ThetaSR%*%t(R1)%*%ginv(PsiSR1)
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

ptm<-(proc.time()-ptm)

Psifinal<-PsiSRfinal%*%t(PsiSRfinal)
print(Psifinal)

indexes <- align_Q( Qin, Qfinal )
# permute
Qfinal <- Qfinal[ , indexes ]
Psifinal <- Psifinal[ indexes, indexes ]
# Q error
err <- rmsd_Q( Qin, Qfinal )
# P error
Ptemp<-Psifinal
Pin2<-Psiin
Ptemp[lower.tri(Ptemp)] <- 0
Pin2[lower.tri(Pin2)] <- 0

message(err,' ',norm((Pin2-Ptemp),"F")*sqrt(2/(K*(K+1))),' ',minf,' ',ptm[3])
