pkgname <- "jmec"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('jmec')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("ecdat")
### * ecdat

flush(stderr()); flush(stdout())

### Name: ecdat
### Title: Simulated clustered data with event, dropout and censoring
### Aliases: ecdat
### Keywords: datasets

### ** Examples

#see ?jmec for example



cleanEx()
nameEx("jmec")
### * jmec

flush(stderr()); flush(stdout())

### Name: jmec
### Title: Fit joint model for event and censoring, using EM-algorithm.
### Aliases: jmec
### Keywords: htest models

### ** Examples


#example data
data(ecdat)

#quick example: few iterations, too optimistic convergence criterion, only part of dataset used
fit1 <- jmec(
	formula.T = Surv(endtime, event) ~ z1 + tr,
	formula.C = Surv(endtime, dropout) ~ z1 + tr,
	clustid = "clust",
	data = ecdat[1:100,],
	Q = 5000,
	nr.nmax = 50,
	stoptol = 1e-5,
	maxit = 5,
	tol = 0.05)
summary(fit1)

#full example: more iterations, realistic convergence criterion, whole dataset used
#(remove comments)
#fit2 <- jmec(
#	formula.T = Surv(endtime, event) ~ z1 + tr,
#	formula.C = Surv(endtime, dropout) ~ z1 + tr,
#	clustid = "clust",
#	data = ecdat,
#	Q = 5000,
#	nr.nmax = 50,
#	stoptol = 1e-5,
#	maxit = 25,
#	tol = 0.005)
#	summary(fit2)




### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
