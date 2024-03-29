% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/porcupine.R
\name{calculate_statistics}
\alias{calculate_statistics}
\title{Determine statistical significance of PCA score of a pathway}
\usage{
calculate_statistics(res_pca_pathway, res_pca_rndm_set)
}
\arguments{
\item{res_pca_pathway}{Result table of pca_pathway function for one pathway}

\item{res_pca_random_subset}{Result table of pca_random function for
one pathway.}
}
\value{
Table of statistics (pvalue and effect size) for each pathway.
}
\description{
Compares the observed PCA score for a pathway to a set of PCA scores of
random gene sets of the same size as pathway. Calculates p-value
and effect size.
}
