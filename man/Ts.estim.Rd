% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ts.estim.R
\name{Ts.estim}
\alias{Ts.estim}
\alias{print.Ts.estim}
\title{Two-Step Estimator}
\usage{
Ts.estim(formula, data, random, all.m.1 = FALSE, D = "UN(1)",
  itermax = 2000, tole = 1e-06)

\method{print}{Ts.estim}(x, ...)
}
\arguments{
\item{formula}{A formula object, with the response on the left of a \code{~} operator, 
and the covariates on the right. The right hand side of the model must also 
include two special terms: a \code{strata} and a \code{cluster} term 
(ex. \code{formula = Y ~ X1 + X2 + X3 + strata(var_strata) + cluster(var_cluster)}).
The \code{strata} and \code{cluster} functions (from the package \pkg{survival}) are
used to identify the stratification and the cluster variables, respectively.}

\item{data}{A data frame (or object coercible by as.data.frame to a data frame) containing the 
variables in the model.}

\item{random}{A formula object, with a blank on the left of a \code{~} operator, 
and, on the right, the covariates with random coefficients among the 
covariate listed in the model \code{formula} (ex. \code{random = ~ X1 + X3}).
The default is to add random coefficients for every covariates listed in 
the model \code{formula}.}

\item{all.m.1}{\code{TRUE} if sum of Y's in all strata is 1, \code{FALSE} otherwise 
(the default). When in doubt use \code{FALSE} (always works, but slower than
necessary if all stratum sums are 1).}

\item{D}{The form of the between-cluster variance-covariance matrix of the regression 
coefficients (matrix D) : either \code{"UN"} for unstructured matrix D or \code{"UN(1)"} 
(the default) for diagonal matrix D.}

\item{itermax}{maximal number of EM iterations (default = 2000)}

\item{tole}{maximal distance between successive EM iterations tolerated
before declaring convergence (default = 0.000001)}

\item{x}{An object, produced by the \code{\link{Ts.estim}} function, to print.}

\item{\dots}{Further arguments to be passed to \code{print.default}.}
}
\value{
\item{beta}{ A vector: the regression coefficients. }

\item{se}{ A vector: the regression coefficients' standard errors. }

\item{vcov}{ A matrix: the variance-covariance matrix of the regression coefficients. }

\item{D}{ A matrix: estimate of the between-cluster variance-covariance matrix of the regression coefficients (matrix D).}

\item{r.effect}{ The random effect estimates. }

\item{coxph.warn}{ A list of character string vectors. If the \code{\link{coxph}} 
                           function generates one or more warnings when fitting the Cox
                           model to a cluster, a copy of these warnings are stored in 
                           \code{coxph.warn$Cluster_name} where \code{Cluster_name} is 
                           the identification value for the cluster in the data set. 
                           A \code{NULL} list element means that \code{\link{coxph}} did not 
                           produce any warnings for that cluster. }

\item{Call}{ The function call.}
}
\description{
Function that computes the two-step estimator proposed in Craiu et al.
(2011) and its \code{print} method.
}
\details{
Calls \code{\link{coxph}} from the package \pkg{survival}.
}
\examples{
# Two ways for specifying the same model
# Data: bison
# Model: covariates forest, biomass and pmeadow
# Random effects in front of forest and biomass
# Main diagonal covariance structure for D (the default)
way1 <- Ts.estim(formula = Y ~ forest  + biomass + pmeadow + 
        strata(Strata) + cluster(Cluster), data = bison, 
        random = ~ forest + biomass)
way1
way2 <- Ts.estim(formula = bison[,3] ~ as.matrix(bison[,c(6,8:9)]) + 
        strata(bison[,2]) + cluster(bison[,1]), data = bison, 
        random = ~ as.matrix(bison[,c(6,8)]))
way2

# Unstructured covariance for D
Fit <- Ts.estim(formula = Y ~ forest  + biomass + pmeadow + 
        strata(Strata) + cluster(Cluster), data = bison, 
        random = ~ forest + biomass, D="UN")
Fit
}
\author{
Radu V. Craiu, Thierry Duchesne, Daniel Fortin and Sophie Baillargeon
}
\references{
Craiu, R.V., Duchesne, T., Fortin, D. and Baillargeon, S.
(2011), Conditional Logistic Regression with Longitudinal Follow-up and
Individual-Level Random Coefficients: A Stable and Efficient Two-Step
Estimation Method, \emph{Journal of Computational and Graphical Statistics}. 
\bold{20}(3), 767-784.
}
\seealso{
\code{\link{ddim}}
}
\keyword{models}

