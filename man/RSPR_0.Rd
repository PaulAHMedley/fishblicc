% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLICC_Model.R
\name{RSPR_0}
\alias{RSPR_0}
\title{Calculate the spawning potential with no fishing}
\usage{
RSPR_0(Galpha, Gbeta, Mk, blicc_ld)
}
\arguments{
\item{Galpha}{Alpha parameter for the Gamma probability density function that
governs growth variability.}

\item{Gbeta}{Rate parameter for the Gamma distribution growth variability
(=Galpha/Linf)}

\item{Mk}{Natural mortality divided by the growth rate K}

\item{blicc_ld}{A standard data list suitable for the model (see function
\link{blicc_dat})}
}
\value{
The spawning potential (a double) when F=0.
}
\description{
The model calculates the spawning potential by multiplying the proportion of
numbers in each length bin with no fishing by the mature biomass per recruit.
}
\examples{
RSPR_0(Galpha=100, Gbeta=100/50, Mk=1.5, blicc_ld=gillnet_ld)

}
