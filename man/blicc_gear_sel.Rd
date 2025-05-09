% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blicc_stan.R
\name{blicc_gear_sel}
\alias{blicc_gear_sel}
\title{Returns data list from \link{blicc_dat} with new gear-selectivity links}
\usage{
blicc_gear_sel(blicc_ld, gear_sel = NULL, gear = NULL, model_name = NULL)
}
\arguments{
\item{blicc_ld}{A standard data list suitable for the model (see function
\link{blicc_dat})}

\item{gear_sel}{A list of integer vectors indexing which \code{sel_fun} are used
for which gear. It can include multiple selectivity functions for each
gear.  If not provided, it is assumed each selectivity function is linked
to each gear in order. Optional.}

\item{gear}{The names or integer vector of the gears being linked to the
selectivity functions}

\item{model_name}{A name for the model (species or fishery). Optional.}
}
\value{
The data object blicc_ld with the new selectivity-gear links.
}
\description{
Sets up the data list with the new links between the selectivity functions
and the gears. This allows the gears to be associated with more than one
selectivity function (mixtures), so increasing the flexibility of the gear
selectivity. The new references are replaced in the data object, which is
then returned.
}
