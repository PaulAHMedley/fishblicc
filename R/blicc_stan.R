# Model Fit Functions -----------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


#' Fit a Bayesian length interval catch curve to a length frequency data sample
#'
#' The model assumes constant recruitment, double-sided normal selectivity,
#' a single natural mortality parameter, and Gamma-distributed von Bertalanffy
#' growth to define expected numbers of fish in pre-defined length bins.
#' The model assumes a negative binomial likelihood function to fit to a single
#' length frequency sample. Selectivity, mortality, asymptotic mean length and
#' scale parameters are fitted to generate the spawning potential ratio as a
#' measure of the stock status.
#'
#' @details
#' The fitted model estimates mortality and survival through sequential length
#' intervals. This is then used to derive abundance and catch in numbers within
#' each length bin. The mortality is required to remain constant within each
#' bin, but otherwise can vary arbitrarily. However, in this implementation
#' mortality is constrained by a parametric model with constant natural
#' mortality by length and a double-sided normal selectivity function for
#' fishing mortality.
#'
#' The model fits 8 parameters. In general, there is insufficient support from
#' a single length frequency sample to estimate these with any precision, and
#' therefore informative priors are used where necessary. The default values
#' for priors can follow recommendations based on life-history invariants
#' (see `blicc_dat`). However the user is still required to provide a prior
#' on the mean maximum length (Linf).
#'
#' The fit starts at the posterior mode which is found using the Stan optimizer.
#' Then a default 4 MCMC chains are initiated and run concurrently. Comparing
#' the chains provides the diagnostics on whether the MCMC has converged or not.
#'
#' It is quite possible to get messages of the type
#' "Exception: neg_binomial_2_lpmf: ..." when the fitting starts. This should
#' not be a problem as long as the message does not persist. If it does either
#' there is a problem with the input data or the MCMC settings. The MCMC
#' parameters can be set using additional parameters for the fit.
#'
#' There are additional parameters for sampling from the MCMC. These can be
#' used to attempt to fix problems. Some, such as low ESS, can be fixed by
#' increasing the number of iterations (in this case, raise ntarget above 2000).
#' Other problems may require consultation with the `rstan::sampling()`
#' documentation. The most likely significant problem would be "divergent" draws
#' which are clearly identified in the stanfit object and in summary
#' diagnostics. These invalidate the MCMC because they cannot be guaranteed as
#' representative of the underlying posterior. This problem might be best
#' addressed in the first instance by increasing the adapt_delta parameter
#' above the default 0.8.
#' For example, the parameter "control=list(adapt_delta=0.9)" can be added
#' to the parameter list on calling the function.
#'
#' @export
#' @param  blicc_ld    A standard data list suitable for the model
#' (see function `blicc_dat()`)
#' @param  ntarget     Target draw for the MCMC
#' @param  nwarmup     Warm up iterations for the Stan MCMC
#' @param  nchain      Number of chains to run in parallel
#' @param  ...         Other arguments passed to `rstan::sampling()`.
#' @return An object of class stanfit returned by `rstan::sampling()`
#' @examples
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0),
#'                 Linf=c(35, 2), NK=50)
#' stf <- blicc_fit(ld, ntarget=100, nwarmup=200, nchain=1) #Quick test
#'
blicc_fit <- function(blicc_ld,
                      ntarget = 2000,
                      nwarmup = 1000,
                      nchain = 4,
                      ...) {

options(mc.cores = parallel::detectCores())  # Needed for parallel chains

  #### remove for the package version
  # Test version
  # stan_fn <- here("source", "BLICC.stan")
  # stmod <<- stan_model(stan_fn, model_name = "BLICC")
  #

  # Find the posterior mode to start
  # # Distributed version
  res <-
    rstan::optimizing(
      stanmodels$BLICC,
      # Test version
      # res <-
      #   optimizing(
      #     stmod,
      #
      data = blicc_ld,
      init = blicc_ini(),
      hessian = TRUE,
      as_vector = FALSE,
      verbose = FALSE,
      iter = 10000,
      refresh = 500,
      tol_obj = 1e-12,
      tol_rel_obj = 1e3,
      tol_grad = 1e-8,
      tol_rel_grad = 1e5,
      tol_param = 1e-8
    )

  niter <- nwarmup + ntarget / nchain
  # Distributed version
  stf <- rstan::sampling(
    stanmodels$BLICC,

    # Test version
    # stf <- stan(
    #   fit = stmod,
    #   file = stan_fn,
    #   model_name = "BLICC",
    #
    data = blicc_ld,
    chains = nchain,
    iter = niter,
    warmup = nwarmup,
    init = blicc_mcmc_ini(nchain, res$par),
    verbose = FALSE,
    ...
  )
  return(stf)
}


#' Generates a data list and initial parameter values for the BLICC model
#'
#' Vectors for the lower bound for each length bin and the number of fish in
#' each length bin and combined with prior parameters into a data list as
#' expected by the BLICC model. Constant bin widths are not necessary for
#' the model and not assumed. Zero frequencies must be explicitly in the
#' fq vector.
#'
#' @details
#' This is a convenience function. The fitting and other functions accept a
#' linked list of values providing necessary information for the BLICC model,
#' including the length bin boundaries and length frequencies. Some basic error
#' checking is carried out and error message is given if the inputs are
#' incorrect.
#'
#' The data list produced as output contains data and parameters suitable for
#' the model fit and producing output. These values can be changed in the list,
#' but must have the same format. The majority of information is the length
#' bin, frequency data and prior parameters (normal for Linf and log-normal
#' for all other parameters). However, the list also contains NK, the number
#' of knots (nodes) used in the quadrature integration. The default is 110,
#' which is safe but extends running times. Lower values risk inaccurate
#' integration but result in faster fitting.
#'
#' It is important that the length frequencies include zeros for bins which
#' contained no fish. The length frequency should be bounded by a single zero
#' bin as the first and last bin in the frequency.
#'
#' @export
#' @param LLB  A vector of lower length boundaries for each length frequency
#' bin. Required.
#' @param fq   A vector, the same length as LLB, containing the frequency data.
#' Zeroes must be included. Required.
#' @param Linf A vector of two values: mean and sd of the prior maximum mean
#' length for the stock. Required.
#' @param Mk   Natural mortality divided by the growth rate K
#' (usually around 1.5). Optional.
#' @param a    The length-weight parameter: a*L^b. Not used in model fitting,
#' but for plots etc. Optional.
#' @param b    A number with the length-weight exponent (a*L^b, usually close
#' to 3.0) Optional.
#' @param L50  The length at 50% maturity, often referred to as "first"
#' maturity, and primarily applies to females. Must be less than Linf, usually
#' around 0.66 Linf, which is assumed if it is not provided. Optional.
#' @param L95  The length at 95% maturity. Must be greater than L50 and less
#' than Linf. A small increment is added L50 if it is not provided. Optional.
#' @param Flat An integer indicating whether selectivity is flat-topped
#' (==0) or not (!=0).  Optional.
#' @param NK   Number of nodes for the Gauss Laguerre quadrature rule. 110 is a
#' safe value, but extends the run time of all calculations.  Optional.
#' @return     A list of vectors and numbers structured suitable for use in the
#' Stan model BLICC.stan.
#' @examples
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2))
#'
blicc_dat <-
  function(LLB,
           fq,
           Linf,
           Mk = NA,
           a = 1,
           b = 3,
           L50 = NA,
           L95 = NA,
           Flat = 1L,
           NK = 110L) {
    # check input vectors, are the same length, and LLB must be in sequence and
    # parameters are in range
    if (!(is.vector(LLB) &
          is.vector(fq)) | (length(LLB) != length(fq))) {
      print(
        "Error: please supply two vectors of same length: 1) lower bound of bins, 2) frequency in each bin, including zeroes."
      )
      return()
    }
    if (any(diff(LLB) <= 0)) {
      print("Error: lower bound of bins must be unique and in ascending order.")
      return()
    }
    if (length(Linf) < 2) {
      print("Error: the Linf vector must contain a mean and sd for the prior.")
      return()
    }
    if (Linf[1] <= min(LLB)) {
      print(paste(
        "Error: Linf must be greater than the lowest bin value:",
        as.character(min(LLB))
      ))
      return()
    }
    if (Linf[2] <= 0) {
      print("Error: Linf prior sd must be greater than zero.")
      return()
    }

    if (b <= 2 | b > 4) {
      print("Error: Length-weight exponent (b) must be greater than 2 and less than 4.")
      return()
    }
    if (is.na(L50)) {
      L50 <- 0.66 * Linf[1]
    } else {
      if (L50 <= 0.2 * Linf[1] | L50 >= Linf[1]) {
        print("Error: Length at 50% maturity must be greater than 0.2*Linf and less than Linf.")
        return()
      }
    }
    if (is.na(L95)) {
      # Ls = -log(1/0.95 - 1)/(L95-L50)
      Ls <- -log(1 / 0.95 - 1) / (0.05 * (Linf[1] - L50) + L50 - L50)
    } else {
      if (L95 <= L50 | L95 >= Linf[1]) {
        print("Error: Length at 95% maturity must be greater than L50 and less than Linf.")
        return()
      }
      Ls <- -log(1 / 0.95 - 1) / (L95 - L50)
    }
    if (NK < 50) {
      print(
        "Warning: Having fewer than 50 nodes for the Gauss Laguerre quadarture rule is not advised."
      )
    }
    if (!is.na(Mk)) {
      if (Mk <= 0) {
        print("Error: Mk must be greater than zero.")
        return()
      } else {
        lMk <- c(log(Mk), 0.1)
      }
    } else {
      # from Prince et al. 2015
      lMk <-
        c(log(b * (1 - (L50 / Linf[1])) / (L50 / Linf[1])),  0.1)
    }
    # First guess for maximum selectivity
    Smx = c((LLB[which(max(fq) == fq)])[1], 5)      # For double-sided normal

    # Length bin mid points (LMP) used for plotting etc. Not used in Stan model.
    LN <- length(LLB)
    LMP <- c((LLB[-LN] + LLB[-1]) * 0.5, LLB[LN] + 0.5)

    return(
      list(
        NK = NK,
        # Number of knots used in the Gauss-Laguerre quadrature rule
        oBN = length(fq),
        # Number of length bins
        Len = LLB,
        # Lower length boundaries for each bin
        LMP = LMP,
        # Length bin mid points (for plotting etc.)
        fq = fq,
        # Frequency data
        a = a,
        # Length-weight scale
        b = b,
        # Length-weight exponential
        Ls = Ls,
        # Slope parameter for the logistic maturity ogive
        Lm = L50,
        # Length at 50% maturity
        Flat_top = Flat,
        # 0 for flat topped selectivity, !=0 implies double normal.
        poLinfm = Linf[1],
        # Expected Linf
        poLinfs = Linf[2],
        # sd for the normal Linf, see above
        polGam = log(1 / 0.1 ^ 2),
        # 10% CV: could be 5% (log(1/0.05^2)) to 30% (log(1/0.3^2)).
        # 30% CV makes length very uninformative however.
        polGas = 0.25,
        # log-normal sd hyper-parameter
        polMkm = lMk[1],
        # see above
        polMks = lMk[2],
        # see above
        polFkm = lMk[1],
        # "Fully exploited"
        polFks = 2.0,
        # Very weak prior with hyper-parameter sd for the log normal
        polSmxm = log(Smx[1]),
        # Prior log-normal mean hyper-parameter for the selectivity
        # mode parameter (see above)
        polSmxs = log(Smx[2]),
        # Prior log-normal sd hyper-parameter for the selectivity mode
        # parameter (see above)
        polSs1m = -0.1,
        # high steepness
        polSs1s = 2.0,
        # Very weak prior with hyperparameter sd for the log normal
        polSs2m = -5.0,
        # Log-normal mean for the right hand selecivity slope parameter
        polSs2s = 2.0,
        # (Weak) Prior with hyperparameter sd for the log normal
        polNB_phim = log(100),
        # Prior log-normal mean parameter for level of overdispersion of
        # observation error: Medium level of overdispersion
        polNB_phis = 0.5
      )
    )
  }


#' Default start parameters centred on priors
#'
#' Provides a list of start parameters for the MCMC - defaults to zero values
#'
#' @return A list of initial parameter values for the BLICC Stan model
#' @noRd
#'
blicc_ini <- function() {
  return(list(
    nLinf   = 0,
    nGalpha = 0,
    nlMk    = 0,
    nlFk    = 0,
    nSmx    = 0,
    nSs1    = 0,
    nSs2    = 0,
    nNB_phi = 0
  ))
}


#' Provides start parameters to the number of chains required for the MCMC run
#'
#' Provides a list of start values for the BLICC model MCMC using function
#' `blicc_ini` if no parameters are provided.
#'
#' @param  pchain A number indicating the number of chains. Default is 4
#' @param  par    A list of start parameter values (optional)
#' @return A list of start parameters equal to the number of chains for
#' the BLICC Stan model
#' @noRd
#'
blicc_mcmc_ini <- function(pchain = 4, par = NULL) {
  if (is.null(par)) {
    return(rep(list(blicc_ini()), pchain))
  } else {
    return(rep(list(par), pchain))
  }
}


#' Calculates the mean log-probability density ratio between two sets of draws
#'
#' Simple Bayes factor estimate based on the mean posterior density.
#' This only works where the data and the likelihood model remain between
#' the two fits. The intended use is to compare a flat-top selectivity (stf0)
#' with a domed selectivity (stf1). In this case, the only change is a single
#' parameter fixed at zero in stf0, but estimated in stf1. This function
#' allows direct comparison of the estimated total posterior
#' probability between the two fits to help determine how much better a domed
#' selectivity fits the data.
#'
#' @export
#' @param stf1,stf0 Stan model output for the numerator and denominator
#' @return A single number estimating the relative probability in support for
#' model stf1 compared to stf2
#' @examples
#' \dontrun{
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0),
#'                 Linf=c(35, 2))
#' slim_domed <- blicc_fit(ld)
#' ld$Flat <- 0L  # Make selectivity flat-topped (logistic or knife-edged)
#' slim_flat <- blicc_fit(ld)
#' posterior_density_ratio(slim_domed, slim_flat)
#'}
posterior_density_ratio <- function(stf1, stf0) {
  p1 <-
    extract(stf1,
            pars = c("lp__"))
  p0 <-
    extract(stf0,
            pars = c("lp__"))

  meanp1 <-
    max(p1$lp__) + log(sum(exp(p1$lp__ - max(p1$lp__))) / length(p1$lp__))
  meanp0 <-
    max(p0$lp__) + log(sum(exp(p0$lp__ - max(p0$lp__))) / length(p0$lp__))

  #print("Bayes factor estimate:")
  return(exp(meanp1 - meanp0))
}
