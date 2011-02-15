\name{Ts.estim}
\alias{Ts.estim}
%\title{Craiu, Duchesne, Fortin and Baillargeon's two-step estimator}
\title{Two-step estimator}
\description{
Function that computes the two-step estimator proposed in Craiu et al. (2011).
}
\usage{
Ts.estim(Donnees, Proposed.model = c(1, 2), random = 1:length(Proposed.model),
         All.m.1 = FALSE, D = "UN(1)", itermax = 2000, tole = 1e-06)
}
\arguments{
  \item{Donnees}{Dataset with cluster number in column 1, stratum number in column 2, \eqn{Y} values in column 3, covariate values in columns 4, 5, ..., \eqn{p}+3}
  \item{Proposed.model}{Lists the covariates that should be in the model, for instance \code{c(1,3,4)} if we want \eqn{x_1}{x1}, \eqn{x_3}{x3} and \eqn{x_4}{x4}.}
  \item{random}{Lists the covariates with random coefficients among the covariates listed in Proposed.model. For instance, if \code{Proposed.model=c(1,3,4)} and \code{random=c(1,2)}, then \eqn{x_1}{x1} and \eqn{x_3}{x3} have random coefficients. The default is to consider all coefficients random.}
  \item{All.m.1}{\code{TRUE} if sum of Y's in all strata is 1, \code{FALSE} otherwise. When in doubt use \code{FALSE} (always works, but slower than necessary if all stratum sums are 1).}
  \item{D}{\code{"UN"} for unstructured matrix D, \code{"UN(1)"} for diagonal matrix D}
  \item{itermax}{maximal number of EM iterations}
  \item{tole}{maximal distance between successive EM iterations tolerated before declaring convergence}
}
\details{
Calls \code{coxph()} from the library survival via the hidden function \code{Fit.models}.
}
\value{
  \item{beta}{ A vector : the regression coefficients. }
  \item{se}{ A vector : the regression coefficients' standard errors. }  
  \item{D}{ The D matrix. }
}
\references{
Craiu, R.V., Duchesne, T., Fortin, D. and Baillargeon, S. (2011), Conditional Logistic Regression with Longitudinal Follow-up and Individual-Level Random Coefficients: A Stable and Efficient Two-Step Estimation Method, \emph{Journal of Computational and Graphical Statistics}, to be published. 
}
\author{
Radu V. Craiu, Thierry Duchesne, Daniel Fortin and Sophie Baillargeon
}
\examples{
data(bison)

dyad <- bison$Cluster
stratum <- bison$Strata
Y <- bison$Y
x1 <- bison$biomass
x2 <- bison$water
x3 <- bison$biomass2

bison.data <- cbind(dyad,stratum,Y,x1,x2,x3)

# Conditional logistic regression
# There are 10 observations per stratum, with
# 2 cases and 8 controls in each stratum

# Model 1: all three covariates x1, x2, x3
# Random effects in front of x1 and x3
# Main diagonal covariance structure for D
Fit1 <- Ts.estim(bison.data, Proposed.model=c(1,2,3),
   random=c(1,3), All.m.1=FALSE, D="UN(1)")

Fit1

# Model 2: only covariates x1, x3
# Random effects in front of x1 and x3
# Unstructured covariance structure for D
Fit2 <- Ts.estim(bison.data, Proposed.model=c(1,3),
   random=c(1,2), All.m.1=FALSE, D="UN")

Fit2
}
\keyword{ models }
