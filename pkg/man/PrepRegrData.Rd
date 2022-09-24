\name{PrepRegrData}
\alias{PrepRegrData}

\title{
 Data prepraration for Bayesian PGLS.
}

\description{
 The function prepares the analysis data and extract the variance-covariance matrix from the phylogeny and matches the variables based on the species names.
}

\usage{PrepRegrData (data, phylo = NULL, phyloDir = NULL, formula = NULL, 
                         weights = NULL, varsForPhylo = NULL, exclSps = NULL, 
                         treeType = "Newick", ...)}

\arguments{
	\item{ data }{Data frame containing the species names, the response variables, and the predictors.}

	\item{phylo }{A phylogenetic tree. If \dQuote{\code{NULL}}, then a file path to a tree needs to be specified through argument \dQuote{\code{phyloDir}}.}
	
	\item{phyloDir }{Path to a phylogenetic tree}

	\item{formula }{An object of class \code{\link{formula}}. A symbolic description of the model to be fitted. See \dQuote{\code{details}}.}

  \item{weights }{Character string indicating the column in \dQuote{\code{data}} that includes the weights for the regression. Default is \code{NULL} when no weights should be included.}
  
  \item{varsForPhylo }{Variables that should be used to match to the phylogeny.}
  
	\item{exclSps }{Character vector of species to be excluded from the analysis.}
	
	\item{treeType }{Type of phylogenetic tree, choices are \dQuote{Nexux} or \dQuote{Newick}. Default is \dQuote{Nexux}. Only necessary when reading tree from file through \dQuote{\code{phyloDir}}.}
	
	 \item{\dots }{Additional arguments passed to function \code{\link{readNexus}}, see \code{\link{readNexus}} documentation.}

}

\details{
The data frame \code{data} needs to include a column with the species names called \dQuote{\code{species}}. If this column is not present, or does not have the proper column name, an error is returned.

Argument \dQuote{\code{phylo}} requires a phologenetic tree as imported by functions \code{\link{read.newick}} or \code{\link{readNexus}} from package \code{\link{phytools}}.


}

\value{
	\item{fill up }{FILL UP}
}

\author{Fernando Colchero \email{colchero@imada.sdu.dk}}

\seealso{FILL UP}

\examples{FILL UP}

\keyword{FILL UP}

