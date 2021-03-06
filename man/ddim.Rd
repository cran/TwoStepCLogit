% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddim.R
\name{ddim}
\alias{ddim}
\alias{print.ddim}
\title{Data Dimension Statistics}
\usage{
ddim(formula, data)

\method{print}{ddim}(x, ...)
}
\arguments{
\item{formula}{A formula object, with the response on the left of a \code{~} operator 
and, on the right hand side, a \code{strata} and a \code{cluster} term 
(ex. \code{formula = Y ~ strata(var_strata) + cluster(var_cluster)}).
The \code{strata} and \code{cluster} functions (from the package survival) are
used to identify the stratification and the cluster variables, respectively.}

\item{data}{A data frame (or object coercible by as.data.frame to a data frame) containing the 
variables in the model.}

\item{x}{An object, produced by the \code{\link{ddim}} function, to print.}

\item{\dots}{Further arguments to be passed to \code{print.default}.}
}
\value{
\item{Sc}{ The number of strata in each cluster.}

\item{Ystat}{ A data.frame with \code{n}, the numbers of observations per stratum 
                      (\eqn{n^c_s}{ncs}), and \code{m}, the sum of the responses per stratum 
                      (\eqn{m^c_s}{mcs}).}
}
\description{
Function that computes dimension statistics for a data set with clusters and strata
and its \code{print} method.
}
\examples{
dimstat <- ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison)
dimstat
}
\seealso{
\code{\link{Ts.estim}}
}

