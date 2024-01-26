\name{RunBayesPGLS}
\alias{RunBayesPGLS}
\alias{RunBayesPGLS.default}

\title{
 Bayesian PGLS function.
}

\description{
 This function runs Bayesian Phylogenetic Generalized Least Squares (PGLS).
}

\usage{\method{RunBayesPGLS}{default} (formula, data, weights = NULL, 
              phylo = NULL, varsForPhylo = NULL, 
              estLambda = TRUE, niter = 30000,  
              burnin = 10001, thinning = 10, nsim, 
              ncpus, exclSps = NULL)}

\arguments{
	\item{formula }{An object of class \code{\link{formula}}. A symbolic description of the model to be fitted. See \dQuote{\code{details}}.}

	\item{data }{Either a data frame with the dependent variable and predictors, or an object of class \dQuote{\code{BayesPGLSdata}}. See \dQuote{\code{details}}.}

  \item{weights }{Character string indicating the column in \dQuote{\code{data}} that includes the weights for the regression. Default is \code{NULL} when no weights should be included.}
  
	\item{phylo }{Optional, phylogenetic tree as imported by package \code{\link{phytools}}. If argument \code{data} is not of class \dQuote{\code{BayesPGLSdata}}, then the phylogenetic tree needs to be provided.}

  \item{varsForPhylo }{Variables that should be used to match to the phylogeny.}

	\item{estLambda }{Logical indicating if Pagel's \code{lambda} should be estimated.}

	\item{niter }{The total number of MCMC steps. Default is 30,000 iterations.}

	\item{ burnin }{The number of iterations for the burn in. Default is 10,000 iterations (see \dQuote{\code{details}})}

	\item{thinning }{The number of skipped MCMC steps to minimize serial autocorrelation. Default is 10.}

	\item{nsim }{A numerical value for the number of simulations to be run (necessary to calculate convergence statistics). Default is 4.}

	\item{ncpus }{A numerical value for the number of CPUS available to run parallel simulations. Default is 4.}

	\item{exclSps }{Character vector indicating which species, if any, should be excluded from the analysis. Default is \code{NA}.}

}

\details{
The function runs Bayesian phylogenetic generalized least squares, combining direct sampling for the regression parameters and the residual variance, and Metropolis-Hastings for Pagel's \eqn{\lambda}. 
}

\value{
  \item{coefficients }{A matrix with estimated coefficients (i.e. mean values per parameter on the thinned sequences after burnin), which includes standard errors, upper and lower 95\% credible intervals, the \dQuote{\code{zeroCoverage}} which provides a lower or upper values of the value of the posterior distribution function for 0, and the potential scale reduction factor for convergence (see \code{Convergence} value below). }
  \item{DIC }{Basic deviance information criterion (DIC) calculations to be used for model selection (Spiegelhalter \emph{et al.} 2002, Celeux \emph{et al.} 2006). Small differences between values should only be used as a reference (see comments in Spiegelhalter \emph{et al.} 2002). If all or some of the simulations failed, then the returned value is \dQuote{\code{Not calculated}}.}
	\item{pars }{Matrix with the converged and thinned sequences of parameter estimates. Used for plotting.}
	\item{likepost }{Matrix with converged and thinned likelihood and posterior values. Used for calculating DIC.}
	\item{runs }{Raw outputs from each of the \code{nsim} runs.}
	\item{fitted }{Fitted values for diagnostics.}
	\item{settings }{MCMC settings used for simulations (e.g., niter, burnin, etc.).}
	\item{form }{Character string reproducing the formula used for analysis.}
	\item{estLambda }{Logical indicating whether Pagel's lambda was estimated.}
	\item{corPred }{Matrix of correlations between predictors. Should be used as reference to avoid collinearity between predictors.}
	\item{potInflObs }{Matrix indicating which observations could have large leverages or could be potential outliers.}
	\item{data }{Object of class \dQuote{\code{BayesPGLSdata}} used for the analysis.}
}

\author{Fernando Colchero \email{fernando_colchero@eva.mpg.de}}

\seealso{\code{\link{PrepRegrData}} to prepare the data for the Bayesian PGLS (built-in in function \code{RunBayesPGLS}).}

\examples{FILL UP}

\keyword{FILL UP}

