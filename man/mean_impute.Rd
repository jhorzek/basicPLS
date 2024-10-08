% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/missing_data_treatment.R
\name{mean_impute}
\alias{mean_impute}
\title{mean_impute}
\usage{
mean_impute(data, weights)
}
\arguments{
\item{data}{data set with missings}

\item{weights}{vector with weights for each person in the data set}
}
\value{
data set with imputed missings
}
\description{
Basic mean imputation for missing data
}
\examples{
library(plsR)
satisfaction_with_missings <- plsR::satisfaction
missings <- matrix(sample(c(TRUE, FALSE),
                   nrow(plsR::satisfaction)*ncol(plsR::satisfaction),
                   prob = c(.1, .9),
                   replace = TRUE),
                   nrow = nrow(plsR::satisfaction),
                   ncol = ncol(plsR::satisfaction))
satisfaction_with_missings[missings] <- NA
mean_impute(data = satisfaction_with_missings,
            weights = rep(1, nrow(satisfaction_with_missings)))
}
