\name{PerformEquating}
\alias{PerformEquating}
\title{Equating parameter estimates of a Rasch model across different contexts.}
\description{
 This functions allows to equate parameter estimates across different cultural, linguistic or environmental contexts, providing coefficients to report item and raw score parameters on the same metric for cross-context comparisons. 
}
\usage{
PerformEquating(listcountry, tol, tol2, mincommon, step)
}
\arguments{
  \item{listcountry}{A list of \code{RM.w} objects.}
  \item{tol}{  Tolerance from the second iteration step on. After the first iteration, 
  a difference of \code{discr = 1+tol} is considered to define common and unique items. \code{discr} is decreased iteration by iteration.}
  \item{tol2}{  Initial tolerance. At the first iteration step, items are considered common if their item severity differ of \code{tol2} from the standard, hence, it is advisable that \code{tol2} is defined sufficieltly large (1 or 2 logits). At each iteration, \code{tol2 = min(tol2 - step, discr)}, where \code{discr} depends on the maximum observed parameter difference.}
  \item{mincommon}{Number of minimum common items. In a specific country, if the number of common items is less than \code{mincommon}, the country is not used to compute the standard.}
  \item{step}{Speed of decreasing of the algorithm tolerance. After the first iteration, \code{tol2 = min(tol2 - step, discr)}, where \code{discr} depends on the maximum observed parameter difference.}
}
\details{
  This function allows to define a common reference scale across different context (depending on, for instance, the different environment, country, language the scale has been administered). In this document, such differences are refered to as \emph{country} effects, but they have in fact a broader, context-related, meaning. In each country, \emph{common} and \emph{unique} item are defined based on the logit-difference with the standard, and such standard is defined based only on the common items. Because of the link between the definition of common items and the standard reference scale, the algorithm needs to start with a very permissive tolerance, and gradually decrease this tolerance during the subsequent iterations.
  
  Let us indicate with \eqn{n} the number of countries, \eqn{k} the number of items and \eqn{b_{ij}} the estimated item parameter for country \eqn{i = 1, \ldots, n} and item \eqn{j = 1, \ldots, k}. The vector \emph{\bold{b}_i} includes all the estimated item parameters for country \eqn{i}, while the vector \emph{\bold{b}_j} contains the estimated parameter for item \eqn{j} in all countries. By construction, \emph{\bold{b}_i} are zero-centered (we will indicate this property as \eqn{\mu_i = 0}), but they can have arbitrary standard deviation (\eqn{\sigma_i}). \eqn{\sigma_i} is related to the scale discrimination power and reliability, but it also reflects a country-specific, and therefore not comparable, metric. 
  
 The equating algortithm is breafly described below.
  
  EQUATING ALGORITHM
  
  \itemize{
  \item Set the iteration counter \eqn{it = 1}
  \item Report the country-specific item parameters on the same metric, i.e. having the same mean and standard deviation. Standardized, country-specific, item parameters (indicated, from here on, as country std.sev.) are obtained by dividing them by their standard deviation \eqn{\sigma_i}.
  \item Assuming that all items are common, define the global reference standard (from here on, ref.sev.) as median of severity parameters across countries.
  \item For each country, calculate the absolute difference between the standardized parameter, and the reference standard (\emph{\bold{d}_{i} = (d_i1, \ldots, d_ik)}), and define the set of \emph{unique} and \emph{common} items for country \eqn{i} as \emph{\bold{U}_i} = \{j:  d_{ij} > \code{tol2}\}, and \emph{\bold{C}_i} = \{j:  {d}_{ij} < \code{tol2}\}.
  \item Define \code{discr = 1+tol}
  \item \emph{WHILE} \code{discr > tol}
  
  \enumerate{
  \item The new metric is defined based on common items in each country, leading to the computation of expansion coefficients \eqn{\alpha_i} and \eqn{\beta_i}:
\deqn{
\beta_i = sd(comm.ref.sev.)/sd(comm.country.std.sev.),
}
and
\deqn{
\alpha_i = mean(comm.ref.sev.)- \beta_i * mean(comm.country.std.sev.)}
Here, \emph{comm.} stands for \emph{common} that reference and country severities are considered only for items defined as \emph{common}.
Therefore, for each country, the new metric is defined in terms of mean and standard deviation of the set of common items.
\item The country-specific parameters are re-scaled and re-centered based on the expansion coefficients as follows:
\deqn{
new country std.sev. = \alpha_i + \beta_i * old country std.sev.
}
\item Calculate the new reference standard (new ref.sev.) based on the median of the new country std. sev. (based only on common items). Based on the new reference standard, the new parameter difference (\emph{\bold{d}_i}), and the new set of unique (\emph{\bold{U}_i}) and common (\emph{\bold{C}_i}) items is defined. 
\item Re-define the tolerance as \code{discr} = max(\code{tol}, max(\emph{\bold{d}_ij'))}, where \emph{j' in \bold{C}^_i}.
\item If \emph{\bold{C}_i} has cardinality less than \code{mincommon}, country \emph{i} is eliminated from the set of countries that contributes to the definition of the reference standard.
  \item Set \eqn{it = it + 1}
    }
    \item end \emph{WHILE}
\item  Define the expansion coefficients on the actual contry-specific severities. This final coefficients will be used to equate the country-specific metric to a common metric, allowing comparable results.
\deqn{
\beta_i final = sd(comm.ref.sev.)/sd(comm.country.sev.),
}
and
\deqn{
\alpha_i final = mean(comm.ref.sev.)- \beta_i final * mean(comm.country.sev.)}
It is important to note that \eqn{\beta_i final} and \eqn{\alpha_i final} consider the actual country severities, and not the standardized ones. 
  }
  }

\value{
  A list with the following elements: 
  \tabular{ll}{
	\code{aa} \tab Vector of the \eqn{\alpha_i}, one for each country, at convergence. \cr
  
	\code{bb} \tab Vector of the \eqn{\beta_i}, one for each country, at convergence. \cr
  
	\code{aaa} \tab Vector of the \eqn{\alpha_i} \emph{final}, one for each country, at convergence. \cr 
  
	\code{bbb} \tab  Vector of the \eqn{\beta_i} \emph{final}, one for each country, at convergence. \cr 
	\code{id.comm} \tab \code{TRUE/FALSE} matrix with \eqn{n} rows and \eqn{k} columns. \code{TRUE} indicates that that item in that country is common. \cr 
	\code{b} \tab Matrix with \eqn{n} rows and \eqn{k} columns containing the estimated item severities for each country. \cr
  \code{b.scaled} \tab Matrix with \eqn{n} rows and \eqn{k} columns containing the scaled item severities for each country at convergence. \cr
  \code{b.tot} \tab Reference standard severity vector at convergence. \cr
  \code{diff.mat} \tab Matrix with \eqn{n} rows and \eqn{k} columns containing the differences between the country-specific severities (\code{b.scaled}) and the reference standard severities (\code{b.tot}). \cr
  \code{items} \tab Names of the items.   \cr
  \code{names}\tab Names of the countries.  \cr
  \code{diff} \tab The tolerance (\code{discr}) at convergence.  \cr
  \code{eqset} \tab Set of countries that contribute to the reference standard definition, at convergence. \cr
  \code{number} \tab Number of common items per country, at convergence.  \cr
  
	}
}


\author{Carlo Cafiero \email{carlo.cafiero@fao.org}}
\seealso{
  \code{\link{RM.w}}, \code{\link{prob.assign}}
}
\examples{
\dontrun{
data(data.FAO_country1)
data(data.FAO_country2)
data(data.FAO_country3)
data(data.FAO_country4)
# Questionnaire data and weights for countries 1, 2, 3 and 4
XX.country1 = data.FAO_country1[,1:8]
XX.country2 = data.FAO_country2[,1:8]
XX.country3 = data.FAO_country3[,1:8]
XX.country4 = data.FAO_country4[,1:8]

wt.country1 = data.FAO_country1$wt
wt.country2 = data.FAO_country2$wt
wt.country3 = data.FAO_country3$wt
wt.country4 = data.FAO_country4$wt

# Fit weighted Rasch for each country separately
rr.country1 = RM.w(XX.country1, wt.country1, country = "Country 1")
rr.country2 = RM.w(XX.country2, wt.country2, country = "Country 2")
rr.country3 = RM.w(XX.country3, wt.country3, country = "Country 3")
rr.country4 = RM.w(XX.country4, wt.country4, country = "Country 4")

# Prepare function input
listcountry = list(rr.country1, rr.country2, rr.country3, rr.country4)

# Perform equating
result = PerformEquating(listcountry, tol = 0.35, tol2 = 2, mincommon = 5, step = 1e-9)

# Expansion factors
alphai = result$aaa   
betai = result$bbb

# Compare prevalence curves on a common metric 
s.levels = seq(-5, 5, length.out = 200)
s.levels.eq = sapply(1:length(listcountry), function(i) alphai[i] + betai[i] * s.levels)

plot(s.levels.eq[,1], prob.assign(listcountry[[1]], sthres = s.levels)$sprob, 
type = "l", col = 1, ylim = c(0, 1), xlab = "Latent trait", ylab = "Cumulative prevalence profile")  

for(i in 2:length(listcountry)) {
	lines(s.levels.eq[,i], prob.assign(listcountry[[i]], sthres = s.levels)$sprob, 
  type = "l", col = i) 
}     
legend("topright", lty = 1, col = 1:length(listcountry), legend = result$names, 
bty = "n", x.intersp = 0.5, cex = 0.8)
text(result$b.tot, 0.2, result$items, srt = 90, cex = 0.6)
}
}


