\name{RM.weights-package}
\alias{RM.weights-package}
\alias{RM.weights-package}
\docType{package}
\title{
Rasch model and extensions for survey data, using Conditional Maximum likelihood (CML). 
}
\description{
This package computes the parameter estimates of a Rasch model for dichotomous and polythomous item responses by using weighted CML estimation, including the possibility of using sampling weights. Procedures to calculate prevalence rates of the phenomenon of interest (latent variable) and to equate item parameters across different context are implemented.
}
\details{
\tabular{ll}{
Package: \tab RM.weights \cr
Type: \tab Package \cr
Version: \tab 1.0 \cr
Date: \tab 2014-11-05 \cr
License: \tab  GPL (>= 2) \cr
}

}
\author{
Carlo Cafiero, Mark Nord, Sara Viviani

Maintainer: Sara Viviani <sara.viviani@fao.org>

}

\keyword{ Conditional maximum likelihood, Rasch model, One-parameter logistic model, Sampling weights, Equating, Probabilistic assignment }

\examples{
\dontrun{
data(data.FAO_country1_2013)
# Questionnaire data and weights
XX = data.FAO_country1_2013[,1:8]
wt = data.FAO_country1_2013$wt
# Fit weighted Rasch
rr = RM.w(XX, wt)
pp = prob.assign(rr, sthres = seq(-5, 5, 0.01))
}
}
