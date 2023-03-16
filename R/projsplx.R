## adapted from ALStructure at https://github.com/StoreyLab/alstructure/blob/master/R/factor.R

## Project a vector onto the simplex
##
## @param y Length-`n` vector.
##
## @return Length-`n` vector which is the projection of `y` onto the `n`-dimensional simplex.
##
## @references
## Chen, Y., and X. Ye. 2011. “Projection Onto A Simplex.” ArXiv E-Prints, January.
projsplx <- function(y) {
    m <- length(y)
    bget <- FALSE
    s <- sort(y, decreasing = TRUE)
    tmpsum <- 0

    for (ii in 1:(m - 1)){
        tmpsum <- tmpsum + s[ii]
        tmax <- (tmpsum - 1) / ii
        if (tmax >= s[ii + 1]) {
            bget <- TRUE
            break
        }
    }
    
    if (!bget)
        tmax <- (tmpsum + s[m] - 1)/m
    
    x <- pmax(y - tmax, rep(0, m))
    return( x )
}
