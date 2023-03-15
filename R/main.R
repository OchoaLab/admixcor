library(MASS)

source('Theta_square_root.R')
source('Objective.R')
source('Initialize.R')
source('projsplx.R')
source('admixcor.R')
#source('align_Q.R')

args = commandArgs(trailingOnly=TRUE);
if(length(args) !=2 )
    stop("Rscript AdmixCor kinship.txt(tab delimited file) K(number of populations)\n")

coances=args[1]
K=as.numeric(args[2])

Theta<-as.matrix(read.table(coances, header = FALSE, sep = "\t")) #coancestry matrix
#Qin<-as.matrix(read.table(args[3], header = FALSE, sep = "\t")) #admixture matrix
#Psiin<-as.matrix(read.table(args[4], header = FALSE, sep = "\t")) #Psi matrix

# main run!
obj <- admixcor( Theta, K )


## indexes <- align_Q( Qin, Qfinal )
## # permute
## Qfinal <- Qfinal[ , indexes ]
## Psifinal <- Psifinal[ indexes, indexes ]
## # Q error
## err <- rmsd_Q( Qin, Qfinal )
## # P error
## Ptemp<-Psifinal
## Pin2<-Psiin
## Ptemp[lower.tri(Ptemp)] <- 0
## Pin2[lower.tri(Pin2)] <- 0

## message(err,' ',norm((Pin2-Ptemp),"F")*sqrt(2/(K*(K+1))),' ',minf,' ',ptm[3])
