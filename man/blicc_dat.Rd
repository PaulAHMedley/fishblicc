% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blicc_stan.R
\name{blicc_dat}
\alias{blicc_dat}
\title{Generates a data list and initial parameter values for the BLICC model}
\usage{
blicc_dat(
  model_name = "fishblicc data model",
  LLB,
  fq,
  Linf,
  sel_fun,
  gear_sel = NULL,
  gear_freq = NULL,
  period_freq = NULL,
  freq_names = NULL,
  gear_names = NULL,
  period_names = NULL,
  Catch = 1,
  Mk = NULL,
  ref_length = -1,
  wt_L = NULL,
  ma_L = NULL,
  a = 0.001,
  b = 3,
  L50 = NULL,
  L95 = NULL,
  NK = NULL
)
}
\arguments{
\item{model_name}{A name for the model (species or fishery). Optional.}

\item{LLB}{A vector of lower length boundaries for each length frequency
bin. Required.}

\item{fq}{A list of vectors, each the same length as LLB, containing the
frequency data for each gear. Zeroes must be included. If only one gear,
can be a vector. Required.}

\item{Linf}{A vector of two values: mean and sd of the prior maximum mean
length for the stock. Required.}

\item{sel_fun}{A vector of selectivity function names (character) or index
(1-4), with 1="logistic", 2="normal", 3="ss_normal", 4="ds_normal". Length
must be the same as fq if \code{gear_sel} is not provided. Required.}

\item{gear_sel}{A list of integer vectors indexing which \code{sel_fun} are used
for which gear. It can include multiple selectivity functions for each
gear.  If not provided, it is assumed each selectivity function is linked
to each gear in order. Optional.}

\item{gear_freq}{An integer vector indexing which data frequencies are from
each gear. If not provided, one-to-one link is assumed if possible,
otherwise an error is raised. Optional.}

\item{period_freq}{An integer vector indexing which data frequencies are in
each time period. If not provided, it is assumed there is only one time
period for all frequencies. Optional.}

\item{freq_names}{A vector of names for the data frequencies for reference.
Optional.}

\item{gear_names}{A vector of names for the gears for reference. Optional.}

\item{period_names}{A vector of names for the time periods for reference.
Optional.}

\item{Catch}{A vector of relative catches, one for each gear. Required if the
number of gears is more than one.}

\item{Mk}{Natural mortality divided by the growth rate K (usually around
1.5). Optional.}

\item{ref_length}{The reference length for the inverse length function if it
is used. (default is NA i.e. constant natural mortality)}

\item{wt_L}{Weight (biomass) for each length bin. Optional.}

\item{ma_L}{Mature biomass for each length bin. Optional.}

\item{a}{The length-weight parameter: a*L^b. Not used in model fitting, but
used to calculate \code{wt_L} if that is not provided. for plots etc. Optional.}

\item{b}{The length-weight exponent (a*L^b, usually close to 3.0). Not used
in model fitting, but used to calculate \code{wt_L} if that is not provided.
Optional.}

\item{L50}{The length at 50\% maturity, often referred to as "first" maturity,
and primarily applies to females. Must be less than \code{Linf}, usually around
0.66 Linf, which is assumed if it is not provided. Optional.}

\item{L95}{The length at 95\% maturity. Must be greater than L50 and less
than Linf. A small increment is added to L50 if it is not provided.
Optional.}

\item{NK}{Number of nodes for the Gauss Laguerre quadrature rule. 110 is a
safe value, but extends the run time of all calculations.  If not provided
it is estimated. Optional.}
}
\value{
A list structured in a manner suitable for use in the Stan model
\link{blicc_mpd} and \link{blicc_fit}.
}
\description{
The function requires vectors for the lower bound for each length bin, a list
of the number of fish in each length bin and prior parameters. This
information is compiled into a data list which can be passed onto the BLICC
model. The minimum requirement is the length bins and frequency data, the
asymptotic mean length (\code{Linf}) and the set of one or more selectivity
functions that will be used in the model. If there is more than one gear, you
also have to provide the relative catches in numbers of fish taken by each
gear. All other parameters are optional, although it is strongly recommended
to use estimates of these parameters wherever they are available.
}
\details{
The fitting and other functions accept a linked list of values
providing necessary information for the BLICC model, including the length
bin boundaries and length frequencies. Some basic error checking is carried
out and an error message is given if the inputs are incorrect.

The data list produced as output contains data and parameters suitable for
the BLICC model fit (\link{blicc_fit}; \link{blicc_mpd}) and other producing output
(e.g. \link{plot_prior}). The majority of information is the length bin,
frequency data and prior parameters (normal for \code{Linf} and log-normal for
all other parameters). A warning is issued for any critical default prior
parameters that are applied. No warning is issued for default prior
parameters (such as Fk, selectivity parameters and negative binomial
\code{NBphi}) that should not be very influential to the fit. Prior parameters
can be subsequently changed using the various \code{blip_} functions:
\link{blip_Linf}, \link{blip_Galpha}, \link{blip_Mk}, \link{blip_Fk}, \link{blip_LH} and
\link{blip_NBphi}.

The list also contains NK, the number of knots (nodes) used in the
quadrature integration. The default is 110, which is safe but extends
running times. Lower values risk inaccurate integration but result in
faster fitting.

The length frequencies should include zeros for bins which contained no
fish. It is recommended that the length frequency should be bounded by a
single zero bin as the first and last bin in the frequency.

The relative catches in numbers of fish are required if there is more than
one gear. Where catches are negligible (e.g. a scientific survey), these
can be set to zero so the fishing mortality is set to zero for this gear.

To allow maximum flexibility, prior hyper-parameters for the selectivity
functions (as log-values for means in the log-normal distribution) are
defined in a single vector. These are linked to each function using an
integer matrix with a row for each selectivity function (gear). This
internal structure is a little complicated but allows mixtures. If you have
one selectivity function for each gear (the default), the data
automatically sets up the correct links. If you are using mixtures, you
need to set up the selectivity functions first and link them to each gear
either using this function \code{blicc_dat} or using \link{blicc_selfun} and
\link{blicc_gear_sel} separately. For mixtures, you will need to manually set
the priors for each mixture function using \link{blip_sel} and \link{blip_mix}.
}
\examples{
ld <- blicc_dat(LLB = 25:35,
                fq = list(c(0,1,2,26,72,66,36,24,12,4,0)),
                Linf = c(35, 2),
                sel_fun = 4,
                gear_names = "Gill net")

}
