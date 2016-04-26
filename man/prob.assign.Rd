\name{prob.assign}
\alias{prob.assign}
\title{Probabilistic assignment of cases to classes of severity along the latent trait.}
\description{
   This function assigns cases to classes of severity along the latent trait. This procedure is useful when different (cultural-geographical-linguistic) contexts are compared in terms of the prevalence of some phenomenon. 
}
\usage{
prob.assign(rr=NULL,  rwthres = NULL, sthres = NULL, eps.a = NULL, ext.distr = NULL, 
flex = list(a=NULL, se.a=NULL, d=NULL, XX=NULL, wt=NULL))
}
\arguments{
  \item{rr}{An object of \code{RM.w} type.}
  \item{rwthres}{Thresholds in terms of raw score, corresponding to which thresholds on the latent trait are calculated.}
  \item{sthres}{Thresholds in terms of latent trait. The probability of being beyond these thresholds is caculated.}
  \item{eps.a}{Tolerance for the algorithm that estimates thresholds on the latent trait corresponding to the specified thresholds in terms of raw score (\code{rwthres}). If left unspecified, the tolerance to convergence is \code{0.001}.}
  \item{ext.distr}{An external probability distribution (for instance, derived by a latent class analysis) corresponding to which thresholds on the latent trait can be estimated. This argument can be used only when \code{length(d) = 2}, in the \code{RM.w} object.}
  \item{flex}{Alternative argument to \code{rr}. If \code{rr=NULL}, than \code{flex} needs
  to be specified as a list with the following arguments: \code{a} (person parameters), \code{se.a} (person parameter standard errors),
  \code{d} (assumption on the extreme raw score parameters; if left unspecified leads to estimate the \code{k}-th raw score parameter at location \code{(k-1)+0.5}), \code{XX} (questionnaire data), \code{wt} (sampling weights; if left unspecified leads to sampling weights equal to one for all subjects).}
}
\details{
The probabilistic assignment procedure is particularly useful when comparing results between different populations that might have interpreted some items of a scale differently. The distribution of the raw score is assumed to be a mixture of Gaussian densities, each centred on the raw score parameters and scaled with the corresponding measurement error. The resulting (complementary) cumulative prevalence, weighted by the proportion of individuals in each raw score, is as follows:
\deqn{
P(x)= \sum_{i=0}^N (1-\int_{-Inf}^x f_i(t)dt) p_i = 
\sum_{i=0}^N (1-F_i(x)) p_i,
}
where \eqn{f_i} is the probability density function (and \eqn{F_i} is the cumulative distribution function) of severity levels centred on raw score parameter \eqn{i} and scaled with the corresponding measurement error, \eqn{p_i} is the proportion observed in raw score \eqn{i}, and \eqn{N} is the maximum number of items. The function \eqn{P(x)} can be used to trace a continuous profile of prevalence on the latent trait for the phenomenon of interest. Working on estimated parameters, instead of classifying individuals based on their raw score, it allows to equate different metrics and it facilitates cross-cultural comparisons.
}
\value{
  A list with the following elements: 
	\tabular{ll}{
	\code{sprob} \tab Estimated weighted probability of being beyond thresholds provided in \code{sthres} (\eqn{P(x)}). \cr
	\code{thres} \tab Thresholds on the latent trait calculated corresponding to thresholds in terms of raw-score specified in \code{rwthres}. \cr 
	\code{f} \tab Probability of being beyond the \code{thres} thresholds divided by raw scores. If more than one assumption on the exreme raw scores is made in \code{rr}, it is going to be a list of a number of elements equal to the number of assumptions made. \cr 
	\code{p} \tab Empirical (weighted) distribution beyond the raw scores specified in \code{rwthres}. If only one assumption on the extreme raw scores is made, then 
  \code{colSums(f)} is approximately equal to \code{p}, where the order of approximation depends on 
  \code{eps.a}.\cr 
	\code{extrathres} \tab Thresholds on the latent variable by considering the (optional) distribution provided in \code{extra.distr}. \cr 
  \code{extraf} \tab Similar to the \code{f} value, but considering the (optional) distribution provided in \code{extra.distr}. \cr 
  \code{f_j} \tab Empirical (weighted) distribution across the raw scores. \cr 
	}
}

\author{Sara Viviani \email{sara.viviani@fao.org}, Mark Nord \email{mark.nord@fao.org}}
\seealso{
  \code{\link{RM.w}}, \code{\link{PerformEquating}}
}
\examples{
\dontrun{
data(data.FAO_country1)
# Questionnaire data and weights
XX.country1 = data.FAO_country1[,1:8]
wt.country1 = data.FAO_country1$wt
# Fit weighted Rasch
rr.country1 = RM.w(XX.country1, wt.country1)

# Thresholds on the latent trait corresponding to a given raw score
pp.country1 = prob.assign(rr.country1, rwthres = c(3, 7))

# Table with prevalences and thresholds
tab = cbind("Raw score" = c(3, 7), "Latent trait" =  pp.country1$thres, 
"Prevalence" = colSums(pp.country1$f))
rownames(tab) = c("Thres 1","Thres 2")
tab

# Pre-defined thresholds on the latent trait
sthresh = seq(-3, 3, 0.01)
pp.country1.2 = prob.assign(rr.country1, sthres = sthresh)
# Plot
plot(sthresh,pp.country1.2$sprob, type = "l", xlab = "Thresholds", 
ylab = "Cumulative prevalence", main = "Country1 food insecurity prevalence")
abline(v = pp.country1$thres[1], lty = 2)
abline(v = pp.country1$thres[2], lty = 2)
text(pp.country1$thres[1], colSums(pp.country1$f)[1], 
"Med.or high prev = 65\%",cex = 0.8, pos = 4)
text(pp.country1$thres[2], colSums(pp.country1$f)[2],
"High prev = 29\%", cex = 0.8, pos = 4)

# More than 2 extremes
# Fit the model
rr.country1.d = RM.w(XX.country1, wt.country1, .d = c(0.2, 0.5, 7.5))
# Probabilistic assignment
pp.country1.d = prob.assign(rr.country1.d, sthres = sthresh)
# Plot
plot(sthresh, pp.country1.d$sprob[[1]], type = "l", xlab = "Thresholds", ylab = "Cumulative prevalence")
lines(sthresh, pp.country1.d$sprob[[2]], col = 2)
}

}
