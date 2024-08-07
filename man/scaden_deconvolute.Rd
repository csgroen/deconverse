% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaden.R
\name{scaden_deconvolute}
\alias{scaden_deconvolute}
\title{scaden (Single-cell assisted deconvolutional network) deconvolution of bulk data using an `screference`}
\usage{
scaden_deconvolute(
  bulk_data,
  scref,
  cache_path = "scaden_results",
  retrain = FALSE,
  ...
)
}
\arguments{
\item{bulk_data}{a matrix of genes-by-samples with bulk mixtures}

\item{scref}{an object of class `screference`}

\item{cache_path}{path to cache the intermediate results}

\item{retrain}{if TRUE, model is re-processed with the `bulk_data` reference and retrained}

\item{...}{passed to scaden_scref if retrain = TRUE}

\item{a}{tibble with deconvolution fractions}
}
\description{
scaden (Single-cell assisted deconvolutional network) deconvolution of bulk data using an `screference`
}
\note{
Reference: Menden, Kevin, Mohamed Marouf, Sergio Oller et al., 2020.
“Deep Learning–Based Cell Composition Analysis from Tissue Expression Profiles.”
Science Advances 6 (30): eaba2619. https://doi.org/10.1126/sciadv.aba2619.

See also: https://github.com/theislab/AutoGeneS
}
