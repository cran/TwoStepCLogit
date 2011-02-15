\name{bison}
\alias{bison}
\docType{data}
\title{
Bison data
}
\description{
Bison data
}
\usage{data(bison)}
\format{
  A data frame with 16818 observations on the following 10 variables.
  \describe{
    \item{\code{Cluster}}{a numeric vector: animal ID}
    \item{\code{Strata}}{a numeric vector: stratum ID}
    \item{\code{Y}}{a numeric vector: 1 for visited locations, 0 otherwise}
    \item{\code{water}}{a numeric vector: indicator covariate}
    \item{\code{agric}}{a numeric vector: indicator covariate}
    \item{\code{forest}}{a numeric vector: indicator covariate}
    \item{\code{meadow}}{a numeric vector: indicator covariate}
    \item{\code{biomass}}{a numeric vector: continuous covariate}
    \item{\code{pmeadow}}{a numeric vector: continuous covariate}
    \item{\code{biomass2}}{a numeric vector: continuous covariate}
  }
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Craiu, R.V., Duchesne, T., Fortin, D. and Baillargeon, S. (2011), Conditional Logistic Regression with Longitudinal Follow-up and Individual-Level Random Coefficients: A Stable and Efficient Two-Step Estimation Method, \emph{Journal of Computational and Graphical Statistics}, to be published. 
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
\keyword{datasets}
