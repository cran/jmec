#MAIN EM ITERATION FUNCTION

#EM function (previously called "huang.em")
jmec <- function(formula.T, formula.C, clustid, data, Q = 5000, nr.nmax = 50, stoptol = 1e-5, maxit = 25, tol = 0.005){
	#note: this function uses only relative tolerance, like emControl {mclust}, default value is 1e-3
	#pre-process the formulas: actively suppress the intercept
		formula.T <- update.formula(formula.T, . ~ . - 1)
		formula.C <- update.formula(formula.C, . ~ . - 1)
	#pre-process the data TO DO: sorting?
		tdlist <- jmec.tempdatfun(formula.T = formula.T, formula.C = formula.C, clustid = clustid, data = data, Q = Q)
####	#computing parameter starting values
##		mod1 <- coxme(Surv(endtime, event) ~ t(tdlist$linpred.T) + (1 | clust), data = tdlist$gendat)
##		mod2 <- coxme(Surv(endtime, cens) ~ t(tdlist$linpred.C) + (1 | clust), data = tdlist$gendat)
##		beta.T.start <- as.numeric(fixef(mod1))
##		beta.C.start <- as.numeric(fixef(mod2))
####		theta.start <- as.numeric((ranef(mod1)[[1]] + ranef(mod2)[[1]]/sqrt(ranef(mod2)[[1]]/ranef(mod1)[[1]]))/2)
####		alpha.start <- as.numeric(sqrt(ranef(mod2)[[1]]/ranef(mod1)[[1]]))
##		theta.start <- 1
##		alpha.start <- 1
##		b.start <- (mod1$frail[[1]] + mod2$frail[[1]]/alpha.start)/2
##		E.exp.b.start <- exp((mod1$frail[[1]] + mod2$frail[[1]]/alpha.start)/2) #I know, this is rough, but there is no way to compute E.exp.b from b when the likelihood is not yet known (E.exp.b. is necessary to compute the likelihood in the first place)
##		E.exp.alphab.start <- exp(alpha.start*(mod1$frail[[1]] + mod2$frail[[1]]/alpha.start)/2) #I know, this is rough, but there is no way to compute E.exp.b from b when the likelihood is not yet known (E.exp.b. is necessary to compute the likelihood in the first place)
	#starting values
		beta.T.start = rep(0, length(attr(terms(formula.T), "term.labels")))
		beta.C.start = rep(0, length(attr(terms(formula.C), "term.labels")))
		alpha.start = 1
		theta.start = 1
		b.start = rep(0, nrow(tdlist$gendat))
		E.exp.b.start = 1
		E.exp.alphab.start = 1
	#prepare list with parameter trace, to store estimates produced by each estimation (first slot contains starting values)
		itdat <- lapply(0:maxit, function(i)list(itnum = i, beta.T = beta.T.start, beta.C = beta.C.start, alpha = alpha.start, theta = theta.start, b = b.start, E.exp.b = E.exp.b.start, E.exp.alphab = E.exp.alphab.start, lambda0.T.cum = NA, lambda0.C.cum = NA, qdraws = NA))
	#make dataframe containing relative changes of the estimates (to determine convergence)
		varnames <- c(paste(rep("T", length(attr(terms(formula.T), "term.labels"))), attr(terms(formula.T), "term.labels"), sep = "."), paste(rep("C", length(attr(terms(formula.C), "term.labels"))), attr(terms(formula.C), "term.labels"), sep = "."), "alpha", "theta")
		changedat <- as.data.frame(matrix(ncol = length(varnames), nrow = maxit))
		colnames(changedat) <- varnames
	#make dataframe for convergence information of Newton-Raphson
		nrdat <- data.frame(
			pse1.nstep = rep(as.numeric(NA), maxit),
			pse1.convergence = rep(as.logical(NA), maxit),
			pse23.nstep = rep(as.numeric(NA), maxit),
			pse23.convergence = rep(as.logical(NA), maxit))
	#convergence indicator
		convergence <- FALSE
	#em iterations
		message("Starting EM iterations...");flush.console()
		for(itnum in 1:maxit){
			em.time <- system.time(em <- jmec.fitter(formula.T = formula.T, formula.C = formula.C, beta.T = itdat[[itnum]]$beta.T, beta.C = itdat[[itnum]]$beta.C, alpha = itdat[[itnum]]$alpha, theta = itdat[[itnum]]$theta, b = itdat[[itnum]]$b, E.exp.b = itdat[[itnum]]$E.exp.b, E.exp.alphab = itdat[[itnum]]$E.exp.alphab, Q = Q, nr.nmax = nr.nmax, stoptol = stoptol, tdlist = tdlist))
			message(paste("   E and M step ", itnum, " done in ", as.character(round(as.list(em.time)$elapsed, 3)), " s.", sep = ""));flush.console()
			itdat[[itnum + 1]]$beta.T = em$beta.T;itdat[[itnum + 1]]$beta.C = em$beta.C;itdat[[itnum + 1]]$alpha = em$alpha;itdat[[itnum + 1]]$theta = em$theta;itdat[[itnum + 1]]$b = em$b;itdat[[itnum + 1]]$E.exp.b = em$E.exp.b;itdat[[itnum + 1]]$E.exp.alphab = em$E.exp.alphab #store estimates
			#new: also store cumulative baseline estimates
				itdat[[itnum + 1]]$lambda0.T.cum = em$lambda0.T.cum
				itdat[[itnum + 1]]$lambda0.C.cum = em$lambda0.C.cum
			#store Newton-Raphson convergence information
				nrdat$pse1.nstep[itnum] <- em$parest1.nstep
				nrdat$pse1.convergence[itnum] <- em$parest1.convergence
				nrdat$pse23.nstep[itnum] <- em$parest23.nstep
				nrdat$pse23.convergence[itnum] <- em$parest23.convergence
			#determine convergence
				changedat[itnum,] <- (c(itdat[[itnum + 1]]$beta.T, itdat[[itnum + 1]]$beta.C, itdat[[itnum + 1]]$alpha, itdat[[itnum + 1]]$theta) - c(itdat[[itnum]]$beta.T, itdat[[itnum]]$beta.C, itdat[[itnum]]$alpha, itdat[[itnum]]$theta))/c(itdat[[itnum]]$beta.T, itdat[[itnum]]$beta.C, itdat[[itnum]]$alpha, itdat[[itnum]]$theta)
				if(max(abs(changedat[itnum,])) <= tol) {convergence <- TRUE;message("   Convergence attained.");break}
				if(itnum == maxit & max(abs(changedat[itnum,])) > tol) {message("   Convergence NOT attained!")}
##				print(max(abs(changedat[itnum,])))
		}
	#put parameter trace in convienient dataframe (and summarize the frailties and functions thereof)
	#(do not put the qdraws and baseline hazard in there)
		pardat <- as.data.frame(do.call("rbind", lapply(1:(maxit + 1), function(i)c(itdat[[i]]$itnum, itdat[[i]]$beta.T, itdat[[i]]$beta.C, itdat[[i]]$alpha, itdat[[i]]$theta))))
		colnames(pardat) <- c("em.it", varnames)
	#select iterations upto convergence or maximum number of iterations
		pardat <- pardat[pardat$em.it <= itnum,]
		changedat <- changedat[1:itnum,]
		nrdat <- nrdat[1:itnum,]
	###############
	#preparations for doing inference
		#trim data from unnecessary variables, paste baseline and cumulative baseline hazards, as well as cluster size and other necessary information
			tdlist$gendat <- cbind(
				tdlist$gendat[, c("clust", "endtime", "event", "cens", "ord")],
				lambda0.T = em$lambda0.T,
				lambda0.T.cum = em$lambda0.T.cum,
				lambda0.C = em$lambda0.C,
				lambda0.C.cum = em$lambda0.C.cum,
				E.exp.b = em$E.exp.b,
				E.exp.alphab = em$E.exp.alphab)
				tdlist$gendat$wkl.T <- as.numeric(exp(crossprod(em$beta.T, tdlist$linpred.T))*tdlist$gendat$E.exp.b)
				tdlist$gendat$wkl.C <- as.numeric(exp(crossprod(em$beta.C, tdlist$linpred.C))*tdlist$gendat$E.exp.alphab)
		#sort gendat in descending order, for efficient computation of score equations
			tdlist$gendat <- tdlist$gendat[order(tdlist$gendat$endtime, decreasing = TRUE),]
			tdlist$gendat$ord2 <- 1:nrow(tdlist$gendat)
			#also sort the linear predictors
				tdlist$linpred.T <- tdlist$linpred.T[,tdlist$gendat$ord]
				tdlist$linpred.C <- tdlist$linpred.C[,tdlist$gendat$ord]
		#only compute once elements that only have to be computed once (independent from bdraw), for efficiency
		#LATER: define these earlier
			#na mailtje Huang inverses lambda0 aangepast (alleen voor event of cens tijden, resp.), missing op NA gezet
			tdlist$gendat$lambda0.T.inv <- with(tdlist$gendat, mapply(function(ev, lam)ifelse(ev == 1, 1/lam, NA), event, lambda0.T))
			tdlist$gendat$lambda0.C.inv <- with(tdlist$gendat, mapply(function(cen, lam)ifelse(cen == 1, 1/lam, NA), cens, lambda0.C))
			tdlist$gendat$lambda0.T.inv.2 <- with(tdlist$gendat, mapply(function(ev, lam)ifelse(ev == 1, -1/(lam^2), NA), event, lambda0.T))
			tdlist$gendat$lambda0.C.inv.2 <- with(tdlist$gendat, mapply(function(cen, lam)ifelse(cen == 1, -1/(lam^2), NA), cens, lambda0.C))
			tdlist$gendat$event.sel <- tdlist$gendat$event == 1 #selection variable for event times
			tdlist$gendat$cens.sel <- tdlist$gendat$cens == 1 #selection variable for cens times
			tdlist$gendat$cp.T <- as.numeric(crossprod(em$beta.T, tdlist$linpred.T))
			tdlist$gendat$cp.C <- as.numeric(crossprod(em$beta.C, tdlist$linpred.C))
			tdlist$matlist1 <- do.call("rbind", lapply(1:nrow(tdlist$gendat), function(k)as.numeric(tdlist$linpred.T[,k]%*%t(tdlist$linpred.T[,k])*tdlist$gendat$lambda0.T.cum[k])))
			tdlist$matlist2 <- do.call("rbind", lapply(1:nrow(tdlist$gendat), function(k)as.numeric(tdlist$linpred.C[,k]%*%t(tdlist$linpred.C[,k])*tdlist$gendat$lambda0.C.cum[k])))
			tdlist$matlist3 <- do.call("rbind", lapply(1:nrow(tdlist$gendat), function(k)as.numeric(tdlist$linpred.C[,k]*tdlist$gendat$lambda0.C.cum[k])))
			tdlist$beta.T.length <- length(em$beta.T)
			tdlist$beta.C.length <- length(em$beta.C)
			tdlist$nevent <- sum(tdlist$gendat$event)
			tdlist$ncens <- sum(tdlist$gendat$cens)
			tdlist$firstder.names <- c(rep("beta.T", length(em$beta.T)), rep("lambda0.T", tdlist$nevent), rep("beta.C", length(em$beta.C)), rep("alpha", 1), rep("lambda0.C", tdlist$ncens), rep("theta", 1)) 
			tdlist$firstder.names.length <- length(tdlist$firstder.names)
	#computation of first derivatives (formulas A.3 - A.8), as a function of the frailty draw for each cluster
	#FOR PREVIOUS VERSIONS OF A.4, A.7, A.10, A.11, A.14, A.16 AND A.17, SEE implementation1_ZZZ41.R
	#FOR DISCUSSION ON CORRECT FORM OF A.4, A.7, A.10, A.11, A.14, A.16 AND A.17, SEE "mailtje2.txt" IN FOLDER "Xuelin_Huang"
		#for checks see implementation1_ZZZ39.R
		firstder.fun <- function(i){
			#i <- 2 #for testing 
			bdraws <- sapply(em$qdraws, function(x)x[i])
			tdlist$gendat$bdraw <- rep(bdraws, times = tdlist$clustsize)[tdlist$gendat$ord]
			exp1 <- exp(tdlist$gendat$cp.T + tdlist$gendat$bdraw)
			der.A.3 <- as.numeric(crossprod(
				t(tdlist$linpred.T), tdlist$gendat$event - tdlist$gendat$lambda0.T.cum*exp1))
			der.A.4 <- tdlist$gendat$lambda0.T.inv[tdlist$gendat$event.sel] - cumsum(exp1)[tdlist$gendat$event.sel] #improved after email Huang 8-Nov-2011
			exp2 <- exp(tdlist$gendat$cp.C + em$alpha*tdlist$gendat$bdraw)
			diff1 <- tdlist$gendat$cens - tdlist$gendat$lambda0.C.cum*exp2
			der.A.5 <- as.numeric(crossprod(t(tdlist$linpred.C), diff1))
			der.A.6 <- sum(tdlist$gendat$bdraw*(diff1))
			der.A.7 <- tdlist$gendat$lambda0.C.inv[tdlist$gendat$cens.sel] - cumsum(exp2)[tdlist$gendat$cens.sel] #improved after email Huang 8-Nov-2011
			der.A.8 <- -0.5*(tdlist$nclust/em$theta - sum(bdraws^2)/(em$theta^2))
			firstder <- c(der.A.3, der.A.4, der.A.5, der.A.6, der.A.7, der.A.8)
			firstder.prod <- firstder %*% t(firstder)
			return(firstder.prod)
		}
		#Compute expectation of matrices, using method of provisional means (see also testing at top of this code)
			message("Computing second term of information matrix...");flush.console()
			pb <- txtProgressBar(min = 0, max = Q, style = 3) #create progress bar
			time2 <- system.time({
				inf.expect.2 <- firstder.fun(1);mean.prev <- inf.expect.2;setTxtProgressBar(pb, 1)
				for(i in 2:Q){
					inf.expect.2 <- mean.prev + (firstder.fun(i) - mean.prev)/i
					mean.prev <- inf.expect.2
					setTxtProgressBar(pb, i)
				}})
			close(pb)
			message(paste("   Done in ", as.character(round(as.list(time2)$elapsed, 3)), " s.", sep = ""));flush.console()
	#computation of second derivatives (formulas A.9 - A.18), as a function of the frailty draw for each cluster
		#for checks see implementation1_ZZZ39.R
		secder.fun <- function(i){
			#i <- 2 #for testing
			bdraws <- sapply(em$qdraws, function(x)x[i])
			tdlist$gendat$bdraw <- rep(bdraws, times = tdlist$clustsize)[tdlist$gendat$ord]
			#der.A.9
				tempvec1 <- exp(tdlist$gendat$cp.T + tdlist$gendat$bdraw)
				der.A.9 <- matrix(
					sapply(1:ncol(tdlist$matlist1), function(i)-1*sum(tdlist$matlist1[,i]*tempvec1)),
					tdlist$beta.T.length,
					tdlist$beta.T.length)
			der.A.10 <- tdlist$gendat$lambda0.T.inv.2[tdlist$gendat$event.sel]
			der.A.11 <- -1*do.call("rbind", lapply(1:tdlist$beta.T.length, function(k)cumsum(tdlist$linpred.T[k,]*tempvec1)[tdlist$gendat$event.sel]))
			#der.A.12
				tempvec2 <- exp(tdlist$gendat$cp.C + em$alpha*tdlist$gendat$bdraw)
				der.A.12 <- matrix(
					sapply(1:ncol(tdlist$matlist2), function(i)-1*sum(tdlist$matlist2[,i]*tempvec2)),
					tdlist$beta.C.length,
					tdlist$beta.C.length)
			der.A.13 <- -1*sum((tdlist$gendat$bdraw^2)*tdlist$gendat$lambda0.C.cum*tempvec2)
			der.A.14 <- tdlist$gendat$lambda0.C.inv.2[tdlist$gendat$cens.sel]
			#der.A.15
				der.A.15 <- matrix(
					sapply(1:ncol(tdlist$matlist3), function(i)-1*sum(tdlist$gendat$bdraw*tdlist$matlist3[,i]*tempvec2)), 1,
					tdlist$beta.C.length)
			der.A.16 <- -1*do.call("rbind", lapply(1:tdlist$beta.C.length, function(k)cumsum(tdlist$linpred.C[k,]*tempvec2)[tdlist$gendat$cens.sel]))
			der.A.17 <- -1*cumsum(tdlist$gendat$bdraw*tempvec2)[tdlist$gendat$cens.sel]
			der.A.18 <- -0.5*(-1*tdlist$nclust/(em$theta^2) + 2*sum(bdraws^2)/(em$theta^3))
			#form complete matrix
				secder <- matrix(data = 0, nrow = tdlist$firstder.names.length, ncol = tdlist$firstder.names.length)
				secder[1:tdlist$beta.T.length, 1:tdlist$beta.T.length] <- der.A.9
				diag(secder)[(tdlist$beta.T.length+1):(tdlist$beta.T.length + tdlist$nevent)] <- der.A.10
				secder[1:tdlist$beta.T.length, (tdlist$beta.T.length+1):(tdlist$beta.T.length + tdlist$nevent)] <- der.A.11
				secder[(tdlist$beta.T.length+1):(tdlist$beta.T.length + tdlist$nevent), 1:tdlist$beta.T.length] <- t(der.A.11)
				secder[with(tdlist, (beta.T.length + nevent + 1):(beta.T.length + nevent + beta.C.length)), with(tdlist, (beta.T.length + nevent + 1):(beta.T.length + nevent + beta.C.length))] <- der.A.12
				secder[with(tdlist, beta.T.length + nevent + beta.C.length + 1), with(tdlist, beta.T.length + nevent + beta.C.length + 1)] <- der.A.13
				diag(secder)[with(tdlist, (beta.T.length + nevent + beta.C.length + 2):(beta.T.length + nevent + beta.C.length + 1 + ncens))] <- der.A.14
				secder[with(tdlist, beta.T.length + nevent + beta.C.length + 1), with(tdlist, (beta.T.length + nevent + 1):(beta.T.length + nevent + beta.C.length))] <- der.A.15
				secder[with(tdlist, (beta.T.length + nevent + 1):(beta.T.length + nevent + beta.C.length)), with(tdlist, beta.T.length + nevent + beta.C.length + 1)] <- t(der.A.15)
				secder[with(tdlist, (beta.T.length + nevent + 1):(beta.T.length + nevent + beta.C.length)), with(tdlist, (beta.T.length + nevent + beta.C.length + 2):(beta.T.length + nevent + beta.C.length + 1 + ncens))] <- der.A.16
				secder[with(tdlist, (beta.T.length + nevent + beta.C.length + 2):(beta.T.length + nevent + beta.C.length + 1 + ncens)), with(tdlist, (beta.T.length + nevent + 1):(beta.T.length + nevent + beta.C.length))] <- t(der.A.16)
				secder[with(tdlist, beta.T.length + nevent + beta.C.length + 1), with(tdlist, (beta.T.length + nevent + beta.C.length + 2):(beta.T.length + nevent + beta.C.length + 1 + ncens))] <- der.A.17
				secder[with(tdlist, (beta.T.length + nevent + beta.C.length + 2):(beta.T.length + nevent + beta.C.length + 1 + ncens)), with(tdlist, beta.T.length + nevent + beta.C.length + 1)] <- t(der.A.17)
				secder[tdlist$firstder.names.length, tdlist$firstder.names.length] <- der.A.18
			return(secder)
		}
		#Compute expectation of matrices, using method of provisional means (see also testing at top of this code)
			message("Computing first term of information matrix...");flush.console()
			pb <- txtProgressBar(min = 0, max = Q, style = 3) #create progress bar
			time1 <- system.time({
				inf.expect.1 <- secder.fun(1);mean.prev <- inf.expect.1;setTxtProgressBar(pb, 1)
				for(i in 2:Q){
					inf.expect.1 <- mean.prev + (secder.fun(i) - mean.prev)/i
					mean.prev <- inf.expect.1
					setTxtProgressBar(pb, i)
				}})
			close(pb)
			message(paste("   Done in ", as.character(round(as.list(time1)$elapsed, 3)), " s.", sep = ""));flush.console()
	#Compute information matrix and covariance matrix, tests etc.
		Imat <- -1*inf.expect.1 - inf.expect.2 #information matrix
		covmat <- solve(Imat)
		#library(corpcor)
		#is.positive.definite(covmat) #(covariance matrix should be positive definite)
		se.est <- sqrt(diag(covmat)[tdlist$firstder.names %in% c("beta.T", "beta.C", "alpha", "theta")])
		qnorm.95.2 <- qnorm(0.975)
		par.est <- pardat[itnum + 1,2:ncol(pardat)]
		par.descr <- 
			data.frame(
				effect = names(par.est),
				coef = as.numeric(t(par.est)),
				HR = exp(as.numeric(t(par.est))),
				se.coef = se.est,
				z = as.numeric(t(par.est))/se.est,
				p = 2*pnorm(abs(as.numeric(t(par.est))/se.est), mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE),
				coef.lower.95 = as.numeric(t(par.est)) - qnorm.95.2*se.est,
				coef.upper.95 = as.numeric(t(par.est)) + qnorm.95.2*se.est,
				HR.lower.95 = exp(as.numeric(t(par.est)) - qnorm.95.2*se.est),
				HR.upper.95 = exp(as.numeric(t(par.est)) + qnorm.95.2*se.est)
			)
		#cbind(par.descr[,1], round(par.descr[,(2:ncol(par.descr))], 3))  #for testing
	#############################
	#output
		res <- list(
			formula.T = formula.T,
			formula.C = formula.C,
			clustid = clustid,
			data = data,
			par.est = par.est,
			par.descr = par.descr,
			par.trace = pardat,
			par.relchange = changedat,
			nrdat = nrdat,
			convergence = convergence,
			it.performed = itnum,
			b = itdat[[itnum + 1]]$b,
			lambda0.T.cum = itdat[[itnum + 1]]$lambda0.T.cum,
			lambda0.C.cum = itdat[[itnum + 1]]$lambda0.C.cum)
			message("Ready.");flush.console()
		class(res) <- "jmec"
		return(res)
	}
