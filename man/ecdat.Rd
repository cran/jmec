\name{ecdat}
\Rdversion{1.1}
\alias{ecdat}
\docType{data}
\title{Simulated clustered data with event, dropout and censoring}

\description{Simulated clustered data with event, dropout and censoring, used as in example in \code{\link{jmec}}.}

\usage{data(ecdat)}
\format{
A data frame with 156 observations on the following 4 variables.
\describe{
\item{\code{id}}{patient ID.}
\item{\code{clust}}{identifier for cluster.}
\item{\code{z1}}{continuous individual-level covariate.}
\item{\code{tr}}{binary exposure.}
\item{\code{endtime}}{individual endtime.}
\item{\code{event}}{binary indicator for event at endtime.}
\item{\code{dropout}}{binary indicator for (informative) dropout at endtime.}
}
}

\details{
These simulated data are used in a detailed in \code{\link{jmec}}. Note that in these data, both "event" and "dropout" are time-to-event variables, as well as censoring (other than dropout). The dropout is informative for the effect of \code{z1} and \code{tr} on the event, but other censoring is non-informative.
}

\references{
Huang X, Wolfe RA. (2002). A Frailty Model for Informative Censoring. Biometrics, 58, 510-520.
}

\author{Willem M. van der Wal \email{w.m.vd.wal@umcutrecht.nl}}

\seealso{\code{\link{jmec}}.}

\examples{
#see ?jmec for example
}
\keyword{datasets}
