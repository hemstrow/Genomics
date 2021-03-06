% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internals.R
\name{apply.snpR.facets}
\alias{apply.snpR.facets}
\title{Apply functions across snpR facets.}
\usage{
apply.snpR.facets(
  x,
  facets = NULL,
  req,
  fun,
  case = "ps",
  par = FALSE,
  ...,
  stats.type = "all",
  response = NULL,
  maf = FALSE,
  interpolate = NULL
)
}
\arguments{
\item{x}{snpRdata object}

\item{facets}{character or NULL, default NULL. Facets to add.}

\item{req}{character. Data type required by fun. See description.}

\item{fun}{function. Function to apply to data.}

\item{case}{character, default "ps". Type of statistic required by fun. See
description.}

\item{par}{numeric or FALSE, default FALSE. Number of parallel computing cores
to use. Works for some cases/reqs.}

\item{...}{other arguments to be passed to fun.}

\item{stats.type}{character, default "all". Other options "pairwise" or
"stats". Type of statistic under to pull under the ps.pf.psf option.}

\item{response}{character, default NULL. Possible response variable name
passed to specific functions.}

\item{maf}{numeric, defualt NULL. Possible maf, passed to some functions.}

\item{interpolate}{character, default NULL. Possible interpolation option,
passed to some functions.}
}
\description{
Internal function to apply functions across snpR data. Facet desicnations
follow rules described in \code{\link{Facets_in_snpR}}. Null or "all" facet
designations follow typical rules.
}
\details{
This function should never be called externally. Raw statistics with metadata
are returned and are not automatically merged into x.

For examples, look at how this function is called in functions such as
calc_pi, calc_pairwise_fst, ect.

Options:

req:

\itemize{ \item{gs: } genotype tables \item{ac: } ac formatted data
\item{meta.gs: } facet, .snp.id metadata cbound genotype tables.
\item{ac.stats: } ac data cbound to stats \item{meta.ac: } ac data cbound to
snp metadata. \item{snpRdata: } subset snpRdata object. }

case:

\itemize{ \item{ps: } stat calculated per snp. \item{ps.pf: } stat calculated
per snp, but split per facet (such as for private alleles, comparisons only
exist within a facet!) \item{facet.pairwise: } stat calculated pairwise
between facets, per snp otherwise. \item{ps.pf.psf: } stat calculated per snp,
but per both sample and snp facet. }
}
\author{
William Hemstrom
}
