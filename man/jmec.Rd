\name{jmec}
\Rdversion{1.1}
\alias{jmec}

\title{Fit joint model for event and censoring, using EM-algorithm.}

\description{Fit joint model for event and censoring, with clustered data, using EM-algorithm, as described in Huang & Wolfe (2002).}

\usage{
jmec(formula.T, formula.C, clustid, data, Q = 5000, nr.nmax = 50,
         stoptol = 1e-5, maxit = 25, tol = 0.005)
}

\arguments{
\item{formula.T}{a \code{Surv()} formula, representing the effect of time-fixed covariates on the event of interest.}
\item{formula.C}{a \code{Surv()} formula, representing the effect of time-fixed covariates on a censoring event of interest, that is possibly informative for the effect of the covariates on the event of interest.}
\item{clustid}{specifies the cluster identifier.}
\item{data}{specifies the dataset containing the outcome variables, covariates and cluster identifier.}
\item{Q}{number of random numbers drawn in the Metropolis-Hastings algorithm.}
\item{nr.nmax}{maximum number of iterations in Newton-Raphson algorithm.}
\item{stoptol}{stopping tolerance in Newton-Raphson algorithm.}
\item{maxit}{maximum number of iterations in EM algorithm.}
\item{tol}{stopping tolerance in EM algorithm.}
}

\details{
See Huang & Wolfe (2002) for technical details. Apart from the (possibly informative) censoring event, observations can be subsequently censored (e.g. administratively), with both indicators for event and censoring equal to \code{0}. Note that the resulting object is of class \code{jmec}, with corresponding \code{summary()} function.}

\value{
\item{formula.T}{formula for event of interest.}
\item{formula.C}{formula for censoring event.}
\item{clustid}{name of cluster identifier.}
\item{data}{data used.}
\item{par.est}{parameter estimates.}
\item{par.descr}{parameter descriptives, including significance tests.}
\item{par.trace}{parameter trace.}
\item{par.relchange}{relative change in parameters.}
\item{nrdat}{Newton-Raphson convergence information.}
\item{convergence}{logical indicator for convergence of EM algorithm.}
\item{it.performed}{number of EM iterations performed.}
\item{b}{estimated cluster-level frailties.}
\item{lambda0.T.cum}{cumulative baseline estimator for event of interest.}
\item{lambda0.C.cum}{cumulative baseline for censoring event}
}

\section{Missing values}{Currently, the exposure variables, covariates and cluster identifier should not contain missing values.}

\references{
Huang X, Wolfe RA. (2002). A Frailty Model for Informative Censoring. Biometrics, 58, 510-520.
}

\author{Willem M. van der Wal \email{w.m.vd.wal@umcutrecht.nl}}

\seealso{\code{\link{ecdat}}}

\examples{

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

}

\keyword{htest}
\keyword{models}
