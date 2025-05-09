% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_Tables.R
\name{blicc_results}
\alias{blicc_results}
\title{Return a tibble containing a summary of the results from a fishblicc fit}
\usage{
blicc_results(blicc_res)
}
\arguments{
\item{blicc_res}{Results from \link{blicc_fit} (stanfit object), \link{blicc_mpd} or
\link{blicc_ref_pts}}
}
\value{
A tibble summarising results of the model fit object.
}
\description{
This function provides a summary of results from the fitted model in a form
suitable for displaying in a table. Note that other packages, such as rstan,
bayesplot and posterior can extract information from the stanfit object and
the \link{blicc_ref_pts} \code{rp_df} 'draws' object in more detailed form.
}
\examples{
blicc_results(trgl_slim)

}
