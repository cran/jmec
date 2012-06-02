#Newton-Raphson implementation
nr.it <- function(inipar, infcn, nmax = 25, stoptol = 1e-05, gradfunc){
	newpar <- inipar
	convergence <- FALSE
	for (i in 0:nmax){
		oldpar <- newpar
		newpar <- oldpar - solve(gradfunc(oldpar), infcn(oldpar))
		if( sqrt(sum((newpar - oldpar)^2)) <= stoptol){convergence <- TRUE;break}
		}
	list(nstep = i, initial = inipar, root = newpar, funcval = infcn(newpar), convergence = convergence)
	}
