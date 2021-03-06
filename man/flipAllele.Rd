% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QC.R
\name{flipAllele}
\alias{flipAllele}
\title{This is the function for flipping alleles}
\usage{
flipAllele(raw.data, raw.data.ori, refaltList, ix.pop, ix.var, log.mat.var,
  correctFlip = TRUE, analyzeRefAltListOnly = TRUE)
}
\arguments{
\item{raw.data}{The input datasets to be considered flipped}

\item{raw.data.ori}{The input datasets to be considered flipped}

\item{refaltList}{The list consists of ref, alt, pos, af and af.diff.max, as well as the option of whether throw away sites with large af.differences checkAF;}

\item{ix.pop}{index of the population}

\item{ix.var}{index of the variant;}

\item{log.mat.var}{The log for QC procedure;}

\item{correctFlip}{Correct for score and covariance matrices for flipped alleles;}
}
\value{
A list consist of modified raw.data, ix.include and log.mat.var
}
\description{
This is the function for flipping alleles
}

