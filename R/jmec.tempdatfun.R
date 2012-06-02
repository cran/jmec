#function for storing temporary data
jmec.tempdatfun <- function(formula.T, formula.C, clustid, data, Q){
	gendat <- data.frame(
		clust = data[,clustid],
		endtime = as.numeric(getResponse(data, formula.T))[1:nrow(data)],
		event = as.numeric(getResponse(data, formula.T))[(nrow(data)+1):(2*nrow(data))],
		cens = as.numeric(getResponse(data, formula.C))[(nrow(data)+1):(2*nrow(data))],
		b.est = as.numeric(NA),
		E.exp.Bk = as.numeric(NA),
		E.exp.alphaBk = as.numeric(NA),
		b.est.update = as.numeric(NA),
		E.exp.Bk.update = as.numeric(NA),
		E.exp.alphaBk.update = as.numeric(NA),
		wkl.T = as.numeric(NA),
		wkl.C = as.numeric(NA),
		lambda0.T = as.numeric(NA),
		lambda0.C = as.numeric(NA),
		lambda0.T.cum = as.numeric(NA),
		lambda0.C.cum = as.numeric(NA))
	linpred.T <- t(model.matrix(formula.T, data = data))
	linpred.C <- t(model.matrix(formula.C, data = data))
	gendat$ord <- 1:nrow(gendat) #record original order of elements
	nclust <- length(unique(gendat[,"clust"])) #determine number of clusters
	return(list(
		gendat = gendat,
		linpred.T = linpred.T,
		linpred.C = linpred.C,
		pct = as.numeric(NA),
		bseq = as.numeric(NA),
		b.v = rep(0, Q),
		dens.v = rep(0, Q),
		b.can.v = rep(0, Q),
		dens.can.v = rep(0, Q),
		a.v = rep(0, Q),
		linpred.T.clust = as.numeric(NA),
		linpred.C.clust = as.numeric(NA),
		lambda0.T.v = as.numeric(NA),
		lambda0.T.cum.v = as.numeric(NA),
		event.v = as.numeric(NA),
		lambda0.C.v = as.numeric(NA),
		lambda0.C.cum.v = as.numeric(NA),
		cens.v = as.numeric(NA),
		mean.dens = as.numeric(NA),
		sd.dens = as.numeric(NA),
		qdraws = rep(list(rep(as.numeric(NA), Q)), nclust),
		qdraws2 = rep(list(rep(as.numeric(NA), Q)), nclust),
		nclust = nclust,
		clustnum = unique(gendat[,"clust"]),
		clustsize = sapply(split(gendat$clust, gendat$clust), function(x)length(x))
		))
}
