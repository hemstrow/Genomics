% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{plot_sfs}
\alias{plot_sfs}
\title{Plot 1 or 2d site frequency spectra.}
\usage{
plot_sfs(sfs, viridis.option = "inferno", log = TRUE)
}
\arguments{
\item{sfs}{matrix or numeric. Either a 2d site frequency spectra stored in a
matrix, with an additional "pops" attribute containing population IDs, such
as c("POP1", "POP2"), where the first pop is the matrix columns and the
second is the matrix rows, or a 1d site frequency spectra stored as a
numeric vector with a similar pops attribute giving the population name.
These objects can be produced from a dadi input file using
\code{\link{make_SFS}}.}

\item{viridis.option}{character, default "inferno". Viridis color scale
option. See \code{\link[ggplot2]{scale_gradient}} for details.}

\item{log}{logical, default TRUE. If TRUE, the number of SNPs in each SFS
cell is log transformed.}
}
\value{
A ggplot2 plot object of the provided SFS.
}
\description{
Plot 1 or 2d site frequency spectra such as those created by
\code{\link{make_SFS}}. Plots are made using ggplot2, and can be freely
modified as is usual for ggplot objects.
}
\details{
The input SFS is either a 2d site frequency spectra stored in a matrix, with
an additional "pops" attribute containing population IDs, such as c("POP1",
"POP2"), where the first pop is the matrix columns and the second is the
matrix rows, or a 1d site frequency spectra stored as a numeric vector with a
similar pops attribute giving the population name. These objects can be
produced from a dadi input file using \code{\link{make_SFS}}.
}
