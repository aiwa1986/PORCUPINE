% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca_scores.R
\name{get_pathway_ind_scores}
\alias{get_pathway_ind_scores}
\title{Extract patient heterogeneity scores on the fist two principal components
for a list of pathways}
\usage{
get_pathway_ind_scores(
  pathways_list,
  reg_net,
  edges,
  scale_data = TRUE,
  center_data = TRUE
)
}
\arguments{
\item{pathways_list}{list of pathways}

\item{reg_net}{Numeric matrix with samples in columns, features in rows}

\item{edges}{Table, containing information on "reg" and "tar"}

\item{scale_data}{Whether to scale the data (TRUE) or not (FALSE).
Default is TRUE.}

\item{center_data}{Logical, whether to center the data (TRUE) or not (FALSE),
Default is TRUE.}
}
\value{
patient heterogeneity scores for each pathway
}
\description{
This function extracts patient heterogeneity scores on the fist
two principal components for a list of pathways
}
