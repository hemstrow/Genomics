% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_functions.R
\name{calc_weighted_stats}
\alias{calc_weighted_stats}
\title{Calculate weighted averages of previously calculated genetic statistics.}
\usage{
calc_weighted_stats(x, facets = NULL, type = "single")
}
\arguments{
\item{x}{snpRdata object.}

\item{facets}{character, default NULL. Facets for which to calculate weighted
stats (broken down by category). See \code{\link{Facets_in_snpR}} for
details.}

\item{type}{character, default "single". Type of statistic to weight:
\itemize{\item{single: } Statistics calculated in a single subfacet, such
as pi. \item{pairwise: } Statistics calculated pairwise between subfacets,
such as Fst. }}
}
\value{
A snpR data object with weighted statistics merged in, accessable via
  \code{\link{get.snpR.stats}} using type = "pop".
}
\description{
Calculates a weighted average for a statistic, weighted by the number of
called genotypes at each locus. Works for single or pairwise statistics (pi,
ho, fst, etc.). Automatically calculates weighted statistic for every
previously calculated statistic.
}
\details{
Weights are calculated using the equation \deqn{ M_{s} = \frac{\sum_{i =
1}^{n} s_{i} * w_{i}}{\sum_{i = 1}^{n} w_{i}}} Where\eqn{n} is the number of
SNPs, \eqn{s_{i}} is the value of the statistic in SNP \eqn{i}, and
\eqn{w_{i}} is the number of times that SNP was genotyped. Note that this
will correct for a range in sequencing depth within samples, but does not
really correct for variance in sequencing depth between populations or other
facet levels.
}
\examples{
# single
x <- calc_basic_snp_stats(stickSNPs, "pop")
x <- calc_weighted_stats(x, "pop")
get.snpR.stats(x, "pop", "pop")

# pairwise
x <- calc_weighted_stats(x, "pop", type = "pairwise") # fst calculated in last step
get.snpR.stats(x, "pop", "pairwise")

}
