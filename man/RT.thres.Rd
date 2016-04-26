\name{RT.thres}
\alias{RT.thres}
\title{Rasch-Thurstone thresholds for trichotomous scales under CML.}
\description{
  Computes the Rasch-Thurstone thresholds for the CML estimation, starting from Rasch thresholds of a trichotomous scale. 
}
\usage{
RT.thres(R.thres)
}
\arguments{
  \item{R.thres}{Vector of rasch thresholds of a trichotomous scale, i.e. points on the latent trait corresponding to which the caracteristic curves of two conscutive categories cross.}
}
\details{
Rasch-Thurstone thresholds are defined as the intersection point, on the latent trait,
between a given category characteristic curve and all the others (cumulated). Formally, these thresholds for a trichotomous scale are comparable with the item severity estimated for a dichotomous scale. Rasch thresholds, usually provided by the Rasch implementation softwares, are defined as intersection points between one category and the consecutive one, and they are not directly comparable to the item severity estimated for a dichotomous scale.
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{b1.rt} \tab R-T threshold for the first category. \cr
	\code{b2.rt} \tab R-T threshold for the second category. \cr 
	\code{b12.rt} \tab Vector of R-T thresholds for the first and second category. \cr 
	}
}

\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}



