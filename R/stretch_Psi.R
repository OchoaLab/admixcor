# shortcuts to these, which have forms with potential for confusion for which parts get transposed
# provided for completeness, though admixcor doesn't actually use it
#' @export
stretch_Psi <- function( Psi, S_inv ) tcrossprod( S_inv %*% Psi, S_inv )
