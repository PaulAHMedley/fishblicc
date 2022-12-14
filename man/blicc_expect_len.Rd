% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_ObsExp.R
\name{blicc_expect_len}
\alias{blicc_expect_len}
\title{Generate a data frame of length-based expected values}
\usage{
blicc_expect_len(blicc_rp, blicc_ld)
}
\arguments{
\item{blicc_rp}{Posterior draws and reference points tibble from
\code{blicc_ref_pts} function.}

\item{blicc_ld}{A data list from the \code{blicc_dat} function.}
}
\value{
A tibble containing fitted values with respect to length
}
\description{
The functions returns a tibble of draws from the posterior MCMC stanfit
object together with a nested tibble for each draw containing the results
(expected values) in terms of length. The BLICC model and parameter values
drawn from an MCMC are used to calculate the expected values for
selectivity, survival, relative population numbers, total mortality and
expected length frequency for each length bin. The resulting table can be
used in various functions to show results. Note that the resulting table
of draws may be large depending on the number of draws.
}
\examples{
lx_df <- blicc_expect_len(eg_rp, eg_ld)
summary(lx_df)

}
