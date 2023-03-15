
projsplx<-function(y)
{
	m = length(y)
	bget = FALSE
	s <- sort(y, decreasing = TRUE)
	tmpsum = 0

	for (ii in 1:(m - 1)){
		tmpsum = tmpsum + s[ii]
	        tmax = (tmpsum - 1) / ii
		if (tmax >= s[ii + 1]){
			bget = TRUE
		        break
		}
	}

	if (!bget){
		tmax = (tmpsum + s[m] - 1)/m
	}

	x <- pmax(y - tmax, rep(0, m));
}
