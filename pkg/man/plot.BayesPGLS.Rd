\name{plot.BayesPGLS}
\alias{plot.BayesPGLS}

\title{
 Plotting BayesPGLS outputs. 
}

\description{
 Draws plots to inspect outputs from function \code{RunBayesPGLS}
}

\usage{plot.BayesPGLS(x, plot.type = "traces")}

\arguments{
	\item{x }{An object of class \dQuote{\code{BayesPGLS}}}

	\item{plot.type }{Character string specifying the type of plot. Options are \dQuote{\code{traces}} (default), \dQuote{\code{density}}, or \dQuote{\code{diagnostics}}}

}

\details{
If \code{plot.type} is set to \dQuote{\code{traces}}, the function plots the parameter traces from the MCMC outputs, which can be used to visually inspect convergence.

If \code{plot.type} is set to \dQuote{\code{density}}, it will draw the posterior densities of the paramaters.

If \code{plot.type} is set to \dQuote{\code{diagnostics}}, it shows residual and qqplots. 
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{RunBayesPGLS}} to run Bayesian PGLS.}

\examples{FILL UP}

\keyword{FILL UP}

