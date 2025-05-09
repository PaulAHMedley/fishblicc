% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_blip.R
\name{blip_sel_auto}
\alias{blip_sel_auto}
\title{Set automatic selectivity prior parameters in data object}
\usage{
blip_sel_auto(blicc_ld, sel_indx = NULL)
}
\arguments{
\item{blicc_ld}{A standard data list suitable for the model (see function
\link{blicc_dat})}

\item{sel_indx}{An integer vector index of the selectivity functions to be
changed. If not provided, all functions will be replaced.}
}
\value{
The data object blicc_ld but with the new selectivity priors
}
\description{
The priors for each function are defined loosely based on the available data
("empirical Bayes"). It is recommended that priors are weakly informative, so
they are set primarily to aid fitting and discourage values outside a
reasonable range.
}
\details{
The prior hyper-parameters are "estimated" if there is a single selectivity
for a gear. The estimates are based on the 50\% (median), 10\% and 90\% quartiles
of the cumulative frequency. This ensures that the selectivity function is
centred on the frequency data and provides a robust estimate of reasonable
priors. Alternatively, setting parameters manually can be done using the
\link{blip_sel} function.

If a selectivity is made up of a mixture of functions, prior parameters are
not estimated and \link{blip_sel} must be used to set the parameters manually. This
is necessary because the user must propose the hypothesis for the mixtures and
it cannot easily be determined from the data what a suitable hypothesis might
be.
}
\examples{
new_ld <- blip_sel_auto(gillnet_ld)

}
