#fitter function
jmec.fitter <- function(formula.T, formula.C, beta.T, beta.C, alpha, theta, b, E.exp.b, E.exp.alphab, Q, nr.nmax, stoptol, tdlist){
	#put previously estimated or starting values of frailties and functions thereof in tdlist
##		#IMPROVEMENT 2012-05-31 PASTDAT ??
##			pastdat <- data.frame(
##				clust = unique(data$clust),
##				b.est = b,
##				E.exp.Bk = E.exp.b,
##				E.exp.alphaBk = E.exp.alphab)
##			tdlist$gendat$b.est <- NULL
##			tdlist$gendat$E.exp.Bk <- NULL
##			tdlist$gendat$E.exp.alphaBk <- NULL
##			tdlist$gendat <- merge(tdlist$gendat, pastdat, by = "clust", all.x = TRUE)
	#put previously estimated or starting values of frailties and functions thereof in tdlist
		tdlist$gendat$b.est <- b
		tdlist$gendat$E.exp.Bk <- E.exp.b
		tdlist$gendat$E.exp.alphaBk <- E.exp.alphab
	#compute Breslow estimates of baseline hazards, as well as cumulative
		tdlist$gendat$wkl.T <- as.numeric(exp(crossprod(beta.T, tdlist$linpred.T))*tdlist$gendat$E.exp.Bk)
		tdlist$gendat$wkl.C <- as.numeric(exp(crossprod(beta.C, tdlist$linpred.C))*tdlist$gendat$E.exp.alphaBk)
		tdlist$gendat <- tdlist$gendat[order(tdlist$gendat$endtime, decreasing = TRUE),] #for efficient computation, sort based on the endtime (xij in the article)
		tdlist$gendat$lambda0.T <- with(tdlist$gendat, event/cumsum(wkl.T))
		tdlist$gendat$lambda0.C <- with(tdlist$gendat, cens/cumsum(wkl.C))
		tdlist$gendat <- tdlist$gendat[order(tdlist$gendat$endtime),] #sort ascending so that cumulative baseline hazards can be computed
		tdlist$gendat$lambda0.T.cum <- cumsum(tdlist$gendat$lambda0.T)
		tdlist$gendat$lambda0.C.cum <- cumsum(tdlist$gendat$lambda0.C)
		tdlist$gendat <- tdlist$gendat[order(tdlist$gendat$ord),] #restore original order
	#Q draws, note that Q is number of draws (argument to fitter function)
		tdlist$pct <- qnorm(0.999, mean = 0, sd = sqrt(theta)) #range used for quick estimation of parameters of candidate distribution for M-H
		tdlist$bseq <- with(tdlist, seq(-pct, pct, pct/100))
		Q.fun <- function(clustnum){
			#store numbers (in tdlist) so that it is not necessary to repeatedly use an index (for speed)
				tdlist$linpred.T.clust <- as.matrix(tdlist$linpred.T[,tdlist$gendat$clust == clustnum])
				tdlist$linpred.C.clust <- as.matrix(tdlist$linpred.C[,tdlist$gendat$clust == clustnum])
				tdlist$lambda0.T.v <- tdlist$gendat$lambda0.T[tdlist$gendat$clust == clustnum]
				tdlist$lambda0.T.cum.v <- tdlist$gendat$lambda0.T.cum[tdlist$gendat$clust == clustnum]
				tdlist$event.v <- tdlist$gendat$event[tdlist$gendat$clust == clustnum]
				tdlist$lambda0.C.v <- tdlist$gendat$lambda0.C[tdlist$gendat$clust == clustnum]
				tdlist$lambda0.C.cum.v <- tdlist$gendat$lambda0.C.cum[tdlist$gendat$clust == clustnum]
				tdlist$cens.v <- tdlist$gendat$cens[tdlist$gendat$clust == clustnum]
			#roughly estimate mean and var of normal candidate distribution (C++ version)
				cdest <- .C("cdest",
					lvT = as.double(c(t(tdlist$linpred.T.clust))),
					lvC = as.double(c(t(tdlist$linpred.C.clust))),
					nobsT = as.integer(ncol(tdlist$linpred.T.clust)),
					nobsC = as.integer(ncol(tdlist$linpred.C.clust)),
					nvarT = as.integer(nrow(tdlist$linpred.T.clust)),
					nvarC = as.integer(nrow(tdlist$linpred.C.clust)),
					cpT = as.double(rep(0, ncol(tdlist$linpred.T.clust))),
					cpC = as.double(rep(0, ncol(tdlist$linpred.C.clust))),
					beta.T = as.double(beta.T),
					beta.C = as.double(beta.C),
					thetaest = as.double(theta),
					alphaest = as.double(alpha),
					bval = as.double(tdlist$bseq),
					bvallength = as.integer(length(tdlist$bseq)),
					lambda0Tv = as.double(tdlist$lambda0.T.v),
					lambda0Cv = as.double(tdlist$lambda0.C.v),
					lambda0Tvcum = as.double(tdlist$lambda0.T.cum.v),
					lambda0Cvcum = as.double(tdlist$lambda0.C.cum.v),
					eventv = as.integer(tdlist$event.v),
					censv = as.integer(tdlist$cens.v),
					l = as.double(rep(0, ncol(tdlist$linpred.T.clust))),
					fb = as.double(rep(0, length(tdlist$bseq))),
					PACKAGE = "jmec"
					)
				tdlist$mean.dens <- sum(cdest$fb*tdlist$bseq)/sum(cdest$fb)
				tdlist$sd.dens <- sqrt(sum(cdest$fb*((tdlist$bseq - tdlist$mean.dens)^2))/sum(cdest$fb))
			#Metropolis-Hastings loop C++
				mh <- .C("mh",
					lvT = as.double(c(t(tdlist$linpred.T.clust))),
					lvC = as.double(c(t(tdlist$linpred.C.clust))),
					nobsT = as.integer(ncol(tdlist$linpred.T.clust)),
					nobsC = as.integer(ncol(tdlist$linpred.C.clust)),
					nvarT = as.integer(nrow(tdlist$linpred.T.clust)),
					nvarC = as.integer(nrow(tdlist$linpred.C.clust)),
					cpT = as.double(rep(0, ncol(tdlist$linpred.T.clust))),
					cpC = as.double(rep(0, ncol(tdlist$linpred.C.clust))),
					beta.T = as.double(beta.T),
					beta.C = as.double(beta.C),
					thetaest = as.double(theta),
					alphaest = as.double(alpha),
					lambda0Tv = as.double(tdlist$lambda0.T.v),
					lambda0Cv = as.double(tdlist$lambda0.C.v),
					lambda0Tvcum = as.double(tdlist$lambda0.T.cum.v),
					lambda0Cvcum = as.double(tdlist$lambda0.C.cum.v),
					eventv = as.integer(tdlist$event.v),
					censv = as.integer(tdlist$cens.v),
					l = as.double(rep(0, ncol(tdlist$linpred.T.clust))),
					bv = as.double(tdlist$b.v),
					densv = as.double(tdlist$dens.v),
					bcanv = as.double(tdlist$b.can.v),
					denscanv = as.double(tdlist$dens.can.v),
					av = as.double(tdlist$a.v),
					rnorm1 = as.double(rnorm(n = 1, mean = tdlist$mean.dens, sd = tdlist$sd.dens)),
					walkv = as.double(rnorm(n = Q, mean = 0, sd = tdlist$sd.dens)),
					runifv = as.double(runif(n = Q)),
					Q = as.integer(Q),
					PACKAGE = "jmec"
					)
			return(mh$bv)
			}
		tdlist$qdraws <- lapply(unique(tdlist$gendat$clust), function(x)Q.fun(clustnum = x))
		tdlist$gendat$b.est.update <- rep(sapply(tdlist$qdraws, function(x)mean(x)), times = tdlist$clustsize)
		tdlist$gendat$E.exp.Bk.update <- rep(sapply(tdlist$qdraws, function(x)mean(exp(x))), times = tdlist$clustsize)
		tdlist$gendat$E.exp.alphaBk.update <- rep(sapply(tdlist$qdraws, function(x)mean(exp(alpha*x))), times = tdlist$clustsize)
	#Sort gendat in descending order, for efficient computation of score equations
		tdlist$gendat <- tdlist$gendat[order(tdlist$gendat$endtime, decreasing = TRUE),]
	#Also sort the linear predictors
		tdlist$linpred.T <- tdlist$linpred.T[,tdlist$gendat$ord]
		tdlist$linpred.C <- tdlist$linpred.C[,tdlist$gendat$ord]
	#Define profile score equation 1, with derivative, and solve it
		pse1 <- function(beta.T.hat){
			wkl.T <- exp(crossprod(beta.T.hat, tdlist$linpred.T))*tdlist$gendat$E.exp.Bk.update
			wkl.T.cumsum <- cumsum(wkl.T)
			zbar.T <- do.call("rbind", lapply(1:length(beta.T), function(i)cumsum(tdlist$linpred.T[i,]*wkl.T)/wkl.T.cumsum))
			return(sapply(1:length(beta.T), function(i)sum(tdlist$gendat$event*(tdlist$linpred.T[i,] - zbar.T[i,]))))
		}
		der.A.9 <- function(beta.T.hat){
			matlist <- lapply(1:nrow(tdlist$gendat), function(i)
				tdlist$linpred.T[,i]%*%t(tdlist$linpred.T[,i])*
				tdlist$gendat$lambda0.T.cum[i]*
				as.numeric(exp(crossprod(beta.T.hat, tdlist$linpred.T[,i]) + tdlist$gendat$b.est.update[i]))
				)
			return(-1*Reduce("+",matlist))
		}
		parest1 <- nr.it(inipar = beta.T, infcn = pse1, nmax = nr.nmax, stoptol = stoptol, gradfunc = der.A.9)
	#Define profile score equations 2 and 3, with derivatives, and solve them together
		tdlist$qdraws2 <- lapply(tdlist$qdraws, function(x)exp(x)) #preparation (for speed)
		pse23 <- function(par.hat){
			beta.C.hat <- par.hat[1:length(beta.C)]
			alpha.hat <- par.hat[length(beta.C) + 1]
			E.exp.alphaBk.hat <- rep(sapply(tdlist$qdraws2, function(x)mean(x^alpha.hat)), times = tdlist$clustsize)[tdlist$gendat$ord]
			wkl.C <- exp(crossprod(beta.C.hat, tdlist$linpred.C))*E.exp.alphaBk.hat
			wkl.C.cumsum <- cumsum(wkl.C)
			zbar.C <- do.call("rbind", lapply(1:length(beta.C), function(i)cumsum(tdlist$linpred.C[i,]*wkl.C)/wkl.C.cumsum))
			E.Bk.exp.alphaBk.hat <- rep(mapply(function(x1, x2)mean(x1*x2^alpha.hat), x1 = tdlist$qdraws, x2 = tdlist$qdraws2), times = tdlist$clustsize)[tdlist$gendat$ord]
			bbar <- cumsum(E.Bk.exp.alphaBk.hat*exp(crossprod(beta.C.hat, tdlist$linpred.C)))/cumsum(wkl.C)
			return(c(sapply(1:length(beta.C), function(i)sum(tdlist$gendat$cens*(tdlist$linpred.C[i,] - zbar.C[i,]))), sum(tdlist$gendat$cens*(tdlist$gendat$b.est.update - bbar))))
		}
		der.A.12 <- function(par.hat){
			beta.C.hat <- par.hat[1:length(beta.C)]
			alpha.hat <- par.hat[length(beta.C) + 1]
			matlist <- lapply(1:nrow(tdlist$gendat), function(i)
				tdlist$linpred.C[,i]%*%t(tdlist$linpred.C[,i])*
				tdlist$gendat$lambda0.C.cum[i]*
				as.numeric(exp(crossprod(beta.C.hat, tdlist$linpred.C[,i]) + alpha.hat*tdlist$gendat$b.est.update[i]))
				)
			return(-1*Reduce("+",matlist))
		}
		der.A.13 <- function(par.hat){
			beta.C.hat <- par.hat[1:length(beta.C)]
			alpha.hat <- par.hat[length(beta.C) + 1]
			matlist <- lapply(1:nrow(tdlist$gendat), function(i)
				tdlist$gendat$b.est.update[i]^2*
				tdlist$gendat$lambda0.C.cum[i]*
				as.numeric(exp(crossprod(beta.C.hat, tdlist$linpred.C[,i]) + alpha.hat*tdlist$gendat$b.est.update[i]))
				)
			return(-1*Reduce("+",matlist))
		}
		der.A.15 <- function(par.hat){
			beta.C.hat <- par.hat[1:length(beta.C)]
			alpha.hat <- par.hat[length(beta.C) + 1]
			matlist <- lapply(1:nrow(tdlist$gendat), function(i)
				tdlist$gendat$b.est.update[i]*tdlist$linpred.C[,i]*
				tdlist$gendat$lambda0.C.cum[i]*
				as.numeric(exp(crossprod(beta.C.hat, tdlist$linpred.C[,i]) + alpha.hat*tdlist$gendat$b.est.update[i]))
				)
			return(-1*Reduce("+",matlist))
		}
		der23 <- function(par.hat){
			hessian23 <- matrix(ncol = length(beta.C) + 1, nrow = length(beta.C) + 1)
			hessian23[1:length(beta.C), 1:length(beta.C)] <- der.A.12(par.hat)
			hessian23[length(beta.C) + 1, length(beta.C) + 1] <- der.A.13(par.hat)
			hessian23[1:length(beta.C), length(beta.C) + 1] <- der.A.15(par.hat)
			hessian23[length(beta.C) + 1, 1:length(beta.C)] <- 0
			return(hessian23)
			} #(combine 3 elements in Hessian for beta.C and alpha)
		parest23 <- nr.it(inipar = c(beta.C, alpha), infcn = pse23, nmax = nr.nmax, stoptol = stoptol, gradfunc = der23)
	#Profile score equation 4: could give problems (due to negative estimate), but it simply reduces to the following (see theta.PDF)
		theta.update <- sum(sapply(tdlist$qdraws, function(x)mean(x^2)))/length(tdlist$qdraws)
	#return data and linear predictors in original order
		tdlist$gendat <- tdlist$gendat[order(tdlist$gendat$ord),]
		tdlist$linpred.T <- tdlist$linpred.T[,order(tdlist$gendat$ord)]
		tdlist$linpred.C <- tdlist$linpred.C[,order(tdlist$gendat$ord)]
	#return list with results (depending on value of qstore)
		res <- list(
			formula.T = formula.T,
			formula.C = formula.C,
			beta.T = parest1$root,
			beta.C = parest23$root[1:length(beta.C)],
			alpha = parest23$root[length(beta.C) + 1],
			theta = theta.update,
			b = tdlist$gendat$b.est.update,
			E.exp.b = tdlist$gendat$E.exp.Bk.update,
			E.exp.alphab = tdlist$gendat$E.exp.alphaBk.update,
			Q = Q,
			stoptol = stoptol,
			parest1.nstep = parest1$nstep,
			parest1.convergence = parest1$convergence,
			parest23.nstep = parest23$nstep,
			parest23.convergence = parest23$convergence,
			lambda0.T = tdlist$gendat$lambda0.T,
			lambda0.C = tdlist$gendat$lambda0.C,
			lambda0.T.cum = tdlist$gendat$lambda0.T.cum,
			lambda0.C.cum = tdlist$gendat$lambda0.C.cum,
			qdraws = tdlist$qdraws
			)
		return(res)
	}
