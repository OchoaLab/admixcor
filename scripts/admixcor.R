# set a seed for reproducible replicates
set.seed(1)

library(genio)
library(tibble)

# load all latest source files
files <- list.files( "../R", pattern = "*.R$", full.names = TRUE )
invisible( sapply( files, source, .GlobalEnv ) )

############
### ARGS ###
############

args <- commandArgs(trailingOnly = TRUE);
if ( length(args) != 8 )
    stop("Rscript admixcor.R <coancestry> <K> <gamma> <delta> <Q_true> <Psi_true> <Q_type> <L_type>\n")

coancestry_file <- args[1]
K <- as.numeric( args[2] )
gamma <- as.numeric( args[3] )
delta <- as.numeric( args[4] )
Q_true_file <- args[5]
Psi_true_file <- args[6]
Q_type <- args[7]
L_type <- args[8]

###########
### RUN ###
###########

# load data to process
Theta <- read_matrix( coancestry_file ) # coancestry matrix
Q_true <- read_matrix( Q_true_file ) # admixture matrix
Psi_true <- read_matrix( Psi_true_file ) # Psi matrix

# main run!
obj <- admixcor( Theta, K, gamma = gamma, delta = delta, Q_type = Q_type, L_type = L_type )
Q_est <- obj$Q
Psi_est <- obj$Psi

# print report!
print( obj$report )

############
### EVAL ###
############

# for final error assessment, align results to true matrices
indexes <- align_Q( Q_true, Q_est )
# permute
Q_est <- Q_est[ , indexes ]
Psi_est <- Psi_est[ indexes, indexes ]
# Q error
err_Q <- rmsd_Q( Q_true, Q_est )
# P error
Ptemp<-Psi_est
Pin2<-Psi_true
Ptemp[lower.tri(Ptemp)] <- 0
Pin2[lower.tri(Pin2)] <- 0
err_Psi <- norm((Pin2-Ptemp),"F")*sqrt(2/(K*(K+1)))

# report errors!
errors <- tibble(
    err_Q = err_Q,
    err_Psi = err_Psi
)
print( errors )
