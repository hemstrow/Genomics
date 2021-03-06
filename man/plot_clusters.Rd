% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_clusters}
\alias{plot_clusters}
\title{PCA, tSNE, and umap plots from snpRdata.}
\usage{
plot_clusters(
  x,
  facets = FALSE,
  plot_type = c("PCA", "tSNE", "umap"),
  check_duplicates = FALSE,
  minimum_percent_coverage = FALSE,
  minimum_genotype_percentage = FALSE,
  interpolation_method = "bernoulli",
  dims = 2,
  initial_dims = 50,
  perplexity = FALSE,
  theta = 0,
  iter = 1000,
  viridis.option = "viridis",
  alt.palette = NULL,
  ncp = NULL,
  ncp.max = 5,
  ...
)
}
\arguments{
\item{x}{snpRdata object.}

\item{facets}{character, default NULL. Categorical sample-level metadata
variables by which to color points. Up to two different sample-specific
facets may be provided. See \code{\link{Facets_in_snpR}} for more details.}

\item{plot_type}{character, default c("PCA", "tSNE", "umap"). Types of plots
to be produced, see description.}

\item{check_duplicates}{logical, default FALSE. Checks for any duplicated
individuals, which will cause errors. Since these rarely exist and
drastically slow down function run-time, this defaults to FALSE.}

\item{minimum_percent_coverage}{numeric, default FALSE. Proportion of samples
a SNP must be sequenced in to be used in generating plots.}

\item{minimum_genotype_percentage}{numeric, default FALSE. Proportion of SNPs
a sample must be sequenced at in order to be used in plots.}

\item{interpolation_method}{character, default "bernoulli". Interpolation
method to use for missing data. Options: \itemize{\item{bernoulli:
}{Interpolated via binomial draw for each allele against minor allele
frequency.} \item{af: }{Interpolated by inserting the expected number of
minor alleles at missing data points given loci minor allele frequencies.}}}

\item{dims}{numeric, default 2. Output dimensionality, default 2.}

\item{initial_dims}{numeric, default 50. The number of dimensions retained in
the initial PCA step during tSNE.}

\item{perplexity}{numeric, default FALSE. Perplexity parameter, by default
found by \code{\link[mmtsne]{hbeta}}, with beta = 1.}

\item{theta}{numeric, default 0. Theta parameter from
\code{\link[Rtsne]{Rtsne}}. Default an exhaustive search.}

\item{iter}{numeric, default 1000. Number of tSNE iterations/umap epochs to
perform.}

\item{viridis.option}{character, default "viridis". Viridis color scale
option to use for significance lines and SNP labels. See
\code{\link[ggplot2]{scale_gradient}} for details.}

\item{alt.palette}{charcter or NULL, default NULL. Optional palette of colors
to use instead of the viridis palette.}

\item{ncp}{numeric or NULL, default NULL. Number of components to consider
for iPCA sn format interpolations of missing data. If null, the optimum
number will be estimated, with the maximum specified by ncp.max. This can
be very slow.}

\item{ncp.max}{numeric, default 5. Maximum number of components to check for
when determining the optimum number of components to use when interpolating
sn data using the iPCA approach.}

\item{...}{Other arguments, passed to \code{\link[Rtsne]{Rtsne}} or
\code{\link[umap]{umap}}.}
}
\value{
A list containing: \itemize{ \item{data: } Raw PCA, tSNE, and umap
  plot data. \item{plots: } ggplot PCA, tSNE, and/or umap plots.} Each of
  these two lists may contain one, two, or three objects, one for each PCA,
  tSNE, or umap plot requested, named "pca" and "tsne", and "umap",
  respectively.
}
\description{
Generate a ggplot cluster plot based on PCA, the Barnes-Hut simulation at
theta>0 implemented in \code{\link[Rtsne]{Rtsne}}, or the Uniform Manifold
Approximation and Projection approach implemented in
\code{\link[umap]{umap}}. Works by conversion to the "sn" format described in
\code{\link{format_snps}} with interpolated missing genotypes.
}
\details{
Cluster plots can be produced via, PCA, tSNE, or umap. The PCA point
coordinates are calculated using \code{\link{prcomp}}. By default, the first
two principal coordinates are plotted. A PC matrix will also be returned for
easy plotting of other PCs. tSNE coordinates are calculated via
\code{\link[Rtsne]{Rtsne}}, which should be consulted to for more details
about this method. Stated simply, tSNE attempts to compress a
multi-dimensional PCA (PCs 1:n) into fewer dimensions while retaining as much
information as possible. As such, a tSNE plot can be seen as a representation
of many different PC axis compressed into a single two-dimensional plot. This
compression process is stochastic, and so plots will vary somewhat between
runs, and multiple runs are recommended. Uniform Manifold Approximation and
Projection (UMAP) coordinates are calculated via \code{\link[umap]{umap}}.
UMAP similarly attempts to reduce multi-dimensional results to a two
dimensional visualization.

For more details on tSNE aruments, \code{\link[Rtsne]{Rtsne}} should be
consulted.

Additional arguments to the UMAP can be also be provided. Additional
information on these arguments can be found in
\code{\link[umap]{umap.defaults}}.

Data points for individuals can be automatically colored by any sample-level
facet categories. Facets should be provided as described in
\code{\link{Facets_in_snpR}}. Up to two different sample-level facets can be
automatically plotted simultaniously.
}
\examples{
# plot colored by population
plot_clusters(stickSNPs, "pop")

# plot colored by population and family
plot_clusters(stickSNPs, "pop", "umap")
}
\references{
Jesse H. Krijthe (2015). Rtsne: T-Distributed Stochastic Neighbor
  Embedding using a Barnes-Hut Implementation, URL:
  \url{https://github.com/jkrijthe/Rtsne}.

Van Der Maaten, L. & Hinton, G. (2008) Visualizing
  high-dimensional data using t-SNE. journal of machine learning research.
  \emph{Journal of Machine Learning Research}.

McInnes, L. & Healy (2018). UMAP: uniform manifold approximation
  and projection. Preprint at URL: \url{https://arxiv.org/abs/1802.03426}.
}
\seealso{
\code{\link[mmtsne]{mmtsne}}
}
\author{
William Hemstrom

Matt Thorstensen
}
