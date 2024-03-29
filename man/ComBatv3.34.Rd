% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ComBatv3.34}
\alias{ComBatv3.34}
\title{ComBat v3.34 batch correction
##' @title ComBat v3.3##' @title ComBat v3.34 batch correction}
\usage{
ComBatv3.34(
  dat,
  batch,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE,
  mean.only = FALSE,
  ref.batch = NULL,
  BPPARAM = BiocParallel::bpparam("SerialParam")
)
}
\arguments{
\item{dat}{Genomic measure matrix (dimensions probe x sample) -
for example, expression matrix}

\item{batch}{Batch covariate (only one batch allowed)}

\item{mod}{Model matrix for outcome of interest and other
covariates besides batch}

\item{par.prior}{(Optional) TRUE indicates parametric adjustments
will be used, FALSE indicates non-parametric adjustments will
be used}

\item{prior.plots}{(Optional) TRUE give prior plots with black as
a kernel estimate of the empirical batch effect density and
red as the parametric}

\item{mean.only}{(Optional) FALSE If TRUE ComBat only corrects the
mean of the batch effect (no scale adjustment)}

\item{ref.batch}{(Optional) NULL If given, will use the selected
batch as a reference for batch adjustment.}

\item{BPPARAM}{(Optional) BiocParallelParam for parallel operation}
}
\value{
data A probe x sample genomic measure matrix, adjusted for batch effects.
}
\description{
Older version of ComBat (sva version 3.34.0). This is for
compatibility with the replication of SCoPE2.

ComBat allows users to adjust for batch effects in datasets where
the batch covariate is known, using methodology described in
Johnson et al. 2007. It uses either parametric or non-parametric
empirical Bayes frameworks for adjusting data for batch effects.
Users are returned an expression matrix that has been corrected
for batch effects. The input data are assumed to be cleaned and
normalized before batch effect removal.
}
