\name{PC.w}
\alias{PC.w}
\title{Estimation of the weighted partial credit model using CML}
\description{
This function computes item and raw score parameter estimates of a partial credit model for both dichotomous and polytomous item responses by using weighted CML estimation. The maximum number of response categories allowed is 4. Input data should be a 0/1 matrix (1 = Yes) for the dichotomous items and 0/1/2... up to the maximum number of categories for the polytomous items. Residual correlation, fit statistics and corresponding standard errors, Rasch reliability and individual fit statistics are also reported.
}
\usage{
PC.w(XX, .wt=NULL, extr=NULL, maxiter=100,minconv=.00001,country=NULL,write.file=F,
            recode = c(0,1,2), write.iteration=F)
}
\arguments{
  \item{XX}{Input 0/1 (for dichotomous items) or 0/1/2/... (for polytomous items) data matrix or data frame; negative responses must be coded as 0s. Rows represent individuals, columns represent items. Missing values are inserted as NA.}
  \item{.wt}{Vector of sampling weights. 
  The length must be the same as the number of rows of \code{XX}. 
  If left unspecified, 
  all the individuals will be equally weighted (\code{.wt = rep(1, nrow(XX)))}. }
  \item{extr}{Optional vector for the assumption on the extreme raw score parameters. Default is \code{0.5} and \code{(k-0.5)},  \code{k} being the maximum number of response categories.} 
    \item{maxiter, minconv}{Convergence criteria.}
    \item{country}{Optional (character) name of the dataset.}
      \item{write.file}{If \code{TRUE}, a CSV file with the main model results will be saved in the working directory.}
      \item{recode}{This is a mandatory argument. It can be 0, 1 or 2. \code{recode = 0} does not imply any category aggregation. \code{recode = 1} aggregates the two more severe categories. \code{recode=2}  aggregates the two less severe categories.}
\item{write.iteration}{If set \code{TRUE} a CSV file with information on the iteration process is produced in the working directory.}
      }
      
\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}
\details{
The weighted CML method is used to estimate the item and Rasch-Thurstone parameters. Respondent parameters  are estimated post-hoc. Cases with missing responses to some items can be included, but will not be used to estimate the Rasch model. 

As the parameters for the extreme raw scores (\code{0} and \code{k}), are undefined under the CML, some assumptions are needed unless the proportions of respondents with those raw scores are so small that they can be considered to be measured without error. 
Vector \code{extr}  gives the possibility to include one assumption for raw score zero and one for raw score \code{k}. The assumption for raw score zero can be a number between 0 and 1 (excluded). The assumption for raw score \code{k} can be a number between \code{k-1} and \code{k} (excluded).
}
\examples{
\dontrun{
data(data.FAO_country1)
# Questionnaire data and weights
XX.country1 = data.FAO_country1[,1:8]
XX.country1[,c(7,8)] = cbind(rbinom(nrow(XX.country1), 3, .1), rbinom(nrow(XX.country1), 3, .2)) 
wt.country1 = data.FAO_country1$wt
# Fit weighted partial credit model with aggregation of the two more severe categories
rr.country1.pol2 = PC.w(XX.country1, wt.country1, recode = 1)
# Write output file
rr.country1.pol2 = PC.w(XX.country1, wt.country1, recode = 1, write.file=T)
}
}

