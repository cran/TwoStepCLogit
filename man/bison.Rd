% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TwoStepCLogit.R
\docType{data}
\name{bison}
\alias{bison}
\title{Bison Dataset}
\format{A data frame with 16818 observations on the following 10 variables.
\describe{ 
  \item{Cluster}{pair of animals (dyad) ID}
  \item{Strata}{stratum ID} 
  \item{Y}{response variable: 1 for visited locations, 0 otherwise} 
  \item{water}{land cover indicator covariate: 1 for water, 0 otherwise} 
  \item{agric}{land cover indicator covariate: 1 for agricultural locations , 0 otherwise} 
  \item{forest}{land cover indicator covariate: 1 for forests, 0 otherwise} 
  \item{meadow}{land cover indicator covariate: 1 for meadows, 0 otherwise}
  \item{biomass}{continuous covariate: above-ground vegetation biomass index measured
                 (in \eqn{kg/m^2}) only at locations within meadows, 0 otherwise}
  \item{pmeadow}{continuous covariate: the proportion of meadow in a circular plot 
                 (700 m in radius) centered at the bison's location}
}}
\description{
Bison data collected in Prince Albert National Park, Saskatchewan, Canada 
(Craiu et al. 2011).
}
\details{
This data set was collected in order to study habitat selection by groups
of free-ranging bison. For each observed group, two individuals (dyad) equipped with GPS 
radio-collars were followed simultaneously. A cluster is defined here as a pair of bison. 
This data set contains 20 clusters. The number of strata per cluster varies
between 13 and 345 for a total of 1410 strata. A stratum is composed of two visited GPS
locations (one for each individual) gathered at the same time, together with 10 
random locations (five drawn within 700 m of each of the two focal bison). Therefore,
there are 12 observations per stratum, with 2 cases (Y=1) and 10 controls (Y=0).
However, due to problems in the data collection, 17 of the 1410 strata have only 6 
observations (1 case and 5 controls).
}
\examples{
# Some descriptive statistics about the data set:
ddim(formula = Y ~ strata(Strata) + cluster(Cluster), data = bison)
 
# Model 1: covariates meadow, biomass and biomass^2
# Random effects in front of biomass and biomass^2
# Main diagonal covariance structure for D
Fit1 <- Ts.estim(formula = Y ~ meadow + biomass + I(biomass^2) + 
        strata(Strata) + cluster(Cluster), data = bison, 
        random = ~ biomass + I(biomass^2), all.m.1=FALSE, D="UN(1)")

Fit1

# Model 2: only covariates biomass and biomass^2
# Random effects in front of biomass and biomass^2
# Main diagonal covariance structure for D
Fit2 <- Ts.estim(formula = Y ~ biomass + I(biomass^2) + strata(Strata) + 
        cluster(Cluster), data = bison, all.m.1=FALSE, D="UN(1)")
Fit2
# Results reported in Table 2 of Craiu et al. (2011).

}
\references{
Craiu, R.V., Duchesne, T., Fortin, D. and Baillargeon, S.
(2011), Conditional Logistic Regression with Longitudinal Follow-up and
Individual-Level Random Coefficients: A Stable and Efficient Two-Step
Estimation Method, \emph{Journal of Computational and Graphical Statistics}. \bold{20}(3), 767-784.
}
\keyword{datasets}

