# BLICC Model Fit Functions -----------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

#' Finds Bayesian length interval catch curve maximum posterior density point
#'
#' The maximum posterior density point is found using on the likelihood for the
#' length frequency data and the priors. The model assumes constant recruitment,
#' double-sided normal selectivity, a single natural mortality parameter, and
#' Gamma-distributed von Bertalanffy growth to define expected numbers of fish
#' in pre-defined length bins. The model also assumes a negative binomial
#' likelihood function to fit to a single length frequency sample. Selectivity,
#' mortality, asymptotic mean length and scale parameters are fitted to generate
#' the spawning potential ratio as a measure of the stock status.
#'
#' It is recommended to test the model and data object using this procedure
#' before investing in the MCMC fit using [blicc_fit].
#'
#' @details Unlike [blicc_fit], the model finds a single maximum posterior point
#'   estimate. This is much faster than the MCMC, but does not provide full
#'   information on uncertainty. The model estimates mortality and survival
#'   through sequential length intervals. This is then used to derive abundance
#'   and catch in numbers within each length bin. The mortality is required to
#'   remain constant within each bin, but otherwise can vary arbitrarily.
#'   However, in this implementation mortality is constrained by a parametric
#'   model with constant natural mortality by length and a double-sided normal
#'   selectivity function for fishing mortality. See [blicc_dat] and [blicc_fit]
#'   documentation for a description of the model and data requirements.
#'
#'   The standard errors in this model are estimated from 1000 random draws from
#'   a multivariate normal with mean at the MPD mode and covariance estimated
#'   from the inverted Hessian matrix. This is an approximation and SE estimates
#'   may differ from the full MCMC sampling.
#'
#'   In addition to the parameters, the fit reports the spawning potential ratio
#'   (SPR) and the log probability at the mode (`lp__`).
#'
#' @export
#' @param  blicc_ld    A standard data list suitable for the model (see function
#'   [blicc_dat])
#' @return A tibble of parameter estimates (mpd) with standard errors.
#' @examples
#' mpd_fit <- blicc_mpd(gillnet_ld)
#' 
blicc_mpd <- function(blicc_ld) {
  # Find the posterior mode
  fit <-
    rstan::optimizing(
      stanmodels$BLICC,
    # Test version
    # rstan::optimizing(
    #   stmod,
      #
      data = blicc_ld,
      init = blicc_ini(blicc_ld),
      hessian = TRUE,
      as_vector = FALSE,
      verbose = FALSE,
      iter = 10000,
      refresh = 500,
      tol_obj = 1e-12,
      tol_rel_obj = 1e3,
      tol_grad = 1e-8,
      tol_rel_grad = 1e5,
      tol_param = 1e-8,
      draws = 1000
    )

  full_par_names <- c("Linf", "Galpha", "Mk",
                      paste0("Fk[", as.character(1:blicc_ld$NF), "]"),
                      paste0("Sm[", as.character(1:(blicc_ld$NP+blicc_ld$NM)), "]"),
                      "NB_phi", "Gbeta", paste0("SPR[", as.character(1:blicc_ld$NT), "]"))
  rd <- fit$theta_tilde[ , full_par_names]

  if (blicc_ld$NF>1) F_names <- paste0("Fk", as.character(1:blicc_ld$NF)) 
  else F_names <- "Fk"
  if (blicc_ld$NT>1) SPR_names <- paste0("SPR", as.character(1:blicc_ld$NT)) 
  else SPR_names <- "SPR"
  
  
  par_names <- c("Linf", "Galpha", "Mk", F_names,
                   paste0("Sm", as.character(1:(blicc_ld$NP+blicc_ld$NM))),
                   "NB_phi", "Gbeta", SPR_names)
  se <- tryCatch(unname(c(apply(rd, MARGIN=2, FUN="sd"), NA)),
                 error=function(cond) return(NA))

  res <- tibble::tibble(
    par = c(full_par_names, "lp__"),
    mpd = unname(c(unlist(fit$par)[par_names],
                   fit$value)),
    se = se
  )
  return(res)
}


#' Fit a Bayesian length interval catch curve to a length frequency data sample
#'
#' The model is fitted to length frequency data using a Markov chain Monte Carlo
#' (MCMC) simulation. The model assumes constant recruitment, a fixed
#' selectivity sub-model, a single natural mortality parameter, and von
#' Bertalanffy growth with gamma probability distribution for individual fish
#' Linf to define expected numbers of fish in pre-defined length bins. The
#' model can fit to multigear data with a length frequency sample from each
#' gear. The model also assumes a negative binomial likelihood function.
#' Selectivity, mortality, asymptotic mean length and scale parameters are
#' fitted to generate the spawning potential ratio as a measure of the stock
#' status.
#'
#' @details The fitted model estimates mortality and survival through sequential
#'   length intervals. This is then used to derive abundance and catch in
#'   numbers within each length bin. The mortality is required to remain
#'   constant within each bin, but otherwise can vary arbitrarily. However, in
#'   this implementation mortality is constrained by parametric models of
#'   selectivity based on the logistic or normal functions with constant natural
#'   mortality by length.
#'
#'   The model fits 4 parameters plus 2-3 parameters for each selectivity
#'   function being fitted. In general, there is insufficient support from a
#'   single length frequency samples to estimate all these with any precision,
#'   and therefore informative priors are used where necessary. The default
#'   values for priors can follow recommendations based on meta-analyses or
#'   life-history invariant (see [blicc_dat]). However the user is still
#'   required to provide a prior on the mean maximum length (Linf) in each case
#'   and this parameter is usually an important pre-requisite for this type of
#'   analysis.
#'
#'   It is recommended to fit a model using the [blicc_mpd] function before
#'   attempting the MCMC fit. This should show up any problems with the model or
#'   data object.
#'
#'   The fit starts at the posterior mode which is found using the Stan
#'   optimizer. Then a default 4 MCMC chains are initiated and run concurrently.
#'   Comparing the chains provides the diagnostics on whether the MCMC has
#'   converged or not.
#'
#'   It is quite possible to get messages of the type "Exception:
#'   neg_binomial_2_lpmf: ..." when the fitting starts. This should not be a
#'   problem as long as the message does not persist. If it does either there is
#'   a problem with the input data or the MCMC settings - particularly the
#'   priors. The priors are setup in the [blicc_dat] function. The MCMC
#'   parameters can be set using additional parameters for the fit.
#'
#'   There are additional parameters for sampling from the MCMC. These can be
#'   used to attempt to fix problems. Some, such as low ESS, can be fixed by
#'   increasing the number of iterations (in this case, raise ntarget above
#'   2000). Other problems may require consultation with the `rstan::sampling()`
#'   documentation. The most likely significant problem would be "divergent"
#'   draws which are clearly identified in the stanfit object and in summary
#'   diagnostics. Significant numbers of these invalidate the MCMC because they
#'   cannot be guaranteed as representative of the underlying posterior.
#'   Divergences are an indication of poor model fit, and it is possible it
#'   indicates this model is not appropriate for the data. The model already
#'   implements non-centred parameterization. Otherwise, it is worth checking
#'   the priors are aligned with the data ([plot_prior]) and then looking for
#'   improvements by adjusting the selectivity model. This problem might also be
#'   addressed by increasing the adapt_delta parameter above the default 0.8,
#'   but this action should probably be a last resort. The parameter
#'   "control=list(adapt_delta=0.9)" can be added to the parameter list on
#'   calling the function.
#'
#'   By default, the start point for the MCMC is the maximum posterior density
#'   (mpd) estimate. These or other estimates can supplied using the `init_fit`
#'   parameter. Start points far from the mpd point may take significantly
#'   longer to converge.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param  ntarget     Target draw for the MCMC
#' @param  nwarmup     Warm up iterations for the Stan MCMC
#' @param  nchain      Number of chains to run in parallel
#' @param  init_fit    List object with element `par` of fitted parameters as
#'   produced from [blicc_mpd] and the `rstan::optimizing` function, for
#'   example. Optional.
#' @param  ...         Other arguments passed to `rstan::sampling()`.
#' @return An object of class stanfit returned by `rstan::sampling()`
#' @examples
#' \dontrun{
#' slim <- blicc_fit(gillnet_ld, ntarget=100, nwarmup=200, nchain=1) #Quick test
#' }
blicc_fit <- function(blicc_ld,
                      ntarget = 2000,
                      nwarmup = 1000,
                      nchain = 4,
                      init_fit = NULL,
                      ...) {

  options(mc.cores = parallel::detectCores())  # Needed for parallel chains

  #### remove for the package version
  # Test version
  # stan_fn <- here("R", "BLICC_mg.stan")
  # stmod <- stan_model(stan_fn, model_name = "BLICC")
  #

  # Find the posterior mode to start
  # # Distributed version
  #if (is.null(init_fit)) {
    init_fit <-
      rstan::optimizing(
        stanmodels$BLICC,
        # Test version
        # rstan::optimizing(
        #   stmod,
        #
        data = blicc_ld,
        init = blicc_ini(blicc_ld),
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
  #}
  
  # New_NK <- LG_Nodes(blicc_ld, rp_df)   # Recalculates the LG knots
  # if (blicc_ld$NK != New_NK) {
  #   glq <- statmod::gauss.quad(New_NK,
  #                              kind = "laguerre", alpha = 0.0)
  #   blicc_ld$NK <- New_NK
  #   blicc_ld$gl_nodes <- glq$nodes
  #   blicc_ld$gl_weights <- glq$weights
  # }
  
  niter <- nwarmup + ntarget / nchain
  # Distributed version
  stf <- rstan::sampling(
    stanmodels$BLICC,

  # Test version
  # stf <- rstan::stan(
  #   fit = stmod,
  #   file = stan_fn,
  #   model_name = "BLICC",
    #
    data = blicc_ld,
    chains = nchain,
    iter = niter,
    warmup = nwarmup,
    init = blicc_mcmc_ini(blicc_ld, nchain, init_fit$par),
    verbose = FALSE,
    ...
  )
  return(stf)
}


#' Generates a data list and initial parameter values for the BLICC model
#'
#' The function requires vectors for the lower bound for each length bin, a list
#' of the number of fish in each length bin and prior parameters. This
#' information is compiled into a data list which can be passed onto the BLICC
#' model. The minimum requirement is the length bins and frequency data, the
#' asymptotic mean length (`Linf`) and the set of one or more selectivity
#' functions that will be used in the model. If there is more than one gear, you
#' also have to provide the relative catches in numbers of fish taken by each
#' gear. All other parameters are optional, although it is strongly recommended
#' to use estimates of these parameters wherever they are available.
#'
#' @details The fitting and other functions accept a linked list of values
#'   providing necessary information for the BLICC model, including the length
#'   bin boundaries and length frequencies. Some basic error checking is carried
#'   out and an error message is given if the inputs are incorrect.
#'
#'   The data list produced as output contains data and parameters suitable for
#'   the BLICC model fit ([blicc_fit]; [blicc_mpd]) and other producing output
#'   (e.g. [plot_prior]). The majority of information is the length bin,
#'   frequency data and prior parameters (normal for `Linf` and log-normal for
#'   all other parameters). A warning is issued for any critical default prior
#'   parameters that are applied. No warning is issued for default prior
#'   parameters (such as Fk, selectivity parameters and negative binomial
#'   `NBphi`) that should not be very influential to the fit. Prior parameters
#'   can be subsequently changed using the various `blip_` functions:
#'   [blip_Linf], [blip_Galpha], [blip_Mk], [blip_Fk], [blip_LH] and
#'   [blip_NBphi].
#'
#'   The list also contains NK, the number of knots (nodes) used in the
#'   quadrature integration. The default is 110, which is safe but extends
#'   running times. Lower values risk inaccurate integration but result in
#'   faster fitting.
#'
#'   The length frequencies should include zeros for bins which contained no
#'   fish. It is recommended that the length frequency should be bounded by a
#'   single zero bin as the first and last bin in the frequency.
#'
#'   The relative catches in numbers of fish are required if there is more than
#'   one gear. Where catches are negligible (e.g. a scientific survey), these
#'   can be set to zero so the fishing mortality is set to zero for this gear.
#'
#'   To allow maximum flexibility, prior hyper-parameters for the selectivity
#'   functions (as log-values for means in the log-normal distribution) are
#'   defined in a single vector. These are linked to each function using an
#'   integer matrix with a row for each selectivity function (gear). This
#'   internal structure is a little complicated but allows mixtures. If you have
#'   one selectivity function for each gear (the default), the data
#'   automatically sets up the correct links. If you are using mixtures, you
#'   need to set up the selectivity functions first and link them to each gear
#'   either using this function `blicc_dat` or using [blicc_selfun] and
#'   [blicc_gear_sel] separately. For mixtures, you will need to manually set
#'   the priors for each mixture function using [blip_sel] and [blip_mix].
#'
#' @export
#' @param model_name A name for the model (species or fishery). Optional.
#' @param LLB  A vector of lower length boundaries for each length frequency
#'   bin. Required.
#' @param fq   A list of vectors, each the same length as LLB, containing the
#'   frequency data for each gear. Zeroes must be included. If only one gear,
#'   can be a vector. Required.
#' @param Linf A vector of two values: mean and sd of the prior maximum mean
#'   length for the stock. Required.
#' @param sel_fun A vector of selectivity function names (character) or index
#'   (1-4), with 1="logistic", 2="normal", 3="ss_normal", 4="ds_normal". Length
#'   must be the same as fq if `gear_sel` is not provided. Required.
#' @param gear_sel A list of integer vectors indexing which `sel_fun` are used
#'   for which gear. It can include multiple selectivity functions for each
#'   gear.  If not provided, it is assumed each selectivity function is linked
#'   to each gear in order. Optional.
#' @param gear_freq An integer vector indexing which data frequencies are from
#'   each gear. If not provided, one-to-one link is assumed if possible,
#'   otherwise an error is raised. Optional.
#' @param period_freq An integer vector indexing which data frequencies are in
#'   each time period. If not provided, it is assumed there is only one time
#'   period for all frequencies. Optional.
#' @param Catch A vector of relative catches, one for each gear. Required if the
#'   number of gears is more than one.
#' @param freq_names A vector of names for the data frequencies for reference.
#'   Optional.
#' @param gear_names A vector of names for the gears for reference. Optional.
#' @param period_names A vector of names for the time periods for reference.
#'   Optional.
#' @param Mk  Natural mortality divided by the growth rate K (usually around
#'   1.5). Optional.
#' @param ref_length  The reference length for the inverse length function if it
#'   is used. (default is NA i.e. constant natural mortality)
#' @param wt_L  Weight (biomass) for each length bin. Optional.
#' @param ma_L  Mature biomass for each length bin. Optional.
#' @param a  The length-weight parameter: a*L^b. Not used in model fitting, but
#'   used to calculate `wt_L` if that is not provided. for plots etc. Optional.
#' @param b  The length-weight exponent (a*L^b, usually close to 3.0). Not used
#'   in model fitting, but used to calculate `wt_L` if that is not provided.
#'   Optional.
#' @param L50 The length at 50% maturity, often referred to as "first" maturity,
#'   and primarily applies to females. Must be less than `Linf`, usually around
#'   0.66 Linf, which is assumed if it is not provided. Optional.
#' @param L95  The length at 95% maturity. Must be greater than L50 and less
#'   than Linf. A small increment is added to L50 if it is not provided.
#'   Optional.
#' @param NK  Number of nodes for the Gauss Laguerre quadrature rule. 110 is a
#'   safe value, but extends the run time of all calculations.  If not provided
#'   it is estimated. Optional.
#' @return  A list structured in a manner suitable for use in the Stan model
#'   [blicc_mpd] and [blicc_fit].
#' @examples
#' ld <- blicc_dat(LLB = 25:35,
#'                 fq = list(c(0,1,2,26,72,66,36,24,12,4,0)),
#'                 Linf = c(35, 2),
#'                 sel_fun = 4,
#'                 gear_names = "Gill net")
#' 
blicc_dat <-
  function(model_name = "fishblicc data model",
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
           NK = NULL) {

    if (!is.vector(LLB, mode="numeric")) {
      stop(
        "Error: LLB required - the lower bound of bins as a vector of unique values and in ascending order. \n"
      )
    }
    if (any(diff(LLB) <= 0)) {
      stop("Error: lower bound of bins (LLB) must be unique and in ascending order. \n")
    }
    NB <- length(LLB)
    LMP <- c((LLB[-NB] + LLB[-1]) * 0.5, LLB[NB] + 0.5*(LLB[NB]-LLB[NB-1])) # Length bin mid points (LMP) used for plotting etc. Not used in Stan model.

    # Data vectors
    if (is.vector(fq, mode="numeric")) fq <- list(fq)  #convert to list if possible
    if (!is.list(fq)) {
      stop(
        "Error: fq (frequency data) must be supplied as a list of vectors (or can be given as a vector if only one.) \n"
      )
    }
    Nfq <- length(fq)
    
    if (Nfq > 1 & all(c(is.null(gear_sel), is.null(gear_freq), is.null(period_freq),
                          is.null(gear_names), is.null(period_names))))
      stop("Error: At least one of gear_sel, gear_freq, period_freq, gear_names or period_names needs to be specified so length frequencies can be allocated correctly./n")    
    
    
    for (i in 1:Nfq)
      if (!is.vector(fq[[i]], mode="numeric") | (NB != length(fq[[i]]))) {
        stop(
          "Error: fq must be a frequency for every bin, including zeroes as vectors of the same length as the number of bins. \n"
        )
      }

    # Parse Time Periods
    ErrStr <- "Error: period_freq must be an integer vector linking each frequency to separate time periods.\n"
    if (is.null(period_freq)) {
      Nperiod <- 1 
      period_freq <- rep(1, Nfq)
    } else { 
      if (length(period_freq) != Nfq)
        stop(ErrStr)
      suppressWarnings(Nperiod <- as.integer(max(period_freq)))
      if(is.na(Nperiod) | Nperiod < 1) 
        stop(ErrStr)
      tp <- sort(unique(period_freq))
      suppressWarnings(
        if (length(tp) != Nperiod | any(tp != 1:Nperiod)) {
          stop(ErrStr)
        })
    }
    if (is.null(period_names)) {
      period_names <- paste0("Period_", as.character(1:Nperiod))
    } else
      if ((length(period_names) != Nperiod) | !is.character(period_names)) {
        stop("Error: period_names must be a character vector with the same length as the number of time periods.\n")
      }
    
    # Parse Gears
    Ngear <- 0
    if (!is.null(gear_freq)) {
      if (length(gear_freq) != Nfq)
        stop(
          "Error: gear_freq must specify the gear applying to each data frequency. \n"
        )
      if (is.character(gear_freq) & is.character(gear_names)) 
        gear_freq <- match(gear_freq, gear_names)
      if (any(is.na(gear_freq))) 
        stop("Error: gear_freq and gear_names must contain identical gear names.\n")
      if (!is.numeric(gear_freq))
          stop("Error: gear_freq must be an vector indexing gears to frequencies.\n")
      Ngear <- max(gear_freq)
      if (length(unique(gear_freq)) != Ngear) 
        stop("Error: gear_freq must use unique identifiers up to the number of gears.\n")
    }
    if (!is.null(gear_sel)) {
      Ngear1 <- length(gear_sel)
      if(is.na(Ngear1) | Ngear1 < 1) 
        stop(
          "Error: gear_sel must be a list of integer vectors that link all selectivity components to each gear.\n"
        )
      if (Ngear > 0 & Ngear1 != Ngear)
        stop("Error: Number of gears referenced in gear_freq and length of gear_sel list must be the same.\n")
      else Ngear <- Ngear1
      if (is.null(gear_freq)) {
        if (Nperiod==1 & Nfq==Ngear) gear_freq <- 1:Nfq 
        else stop("Error: gear_freq needs to be specified indicating which gear produced which length frequency.\n")
      }
    }
    
    if (Ngear==0) {
      if (length(gear_names) == Nfq) {
        Ngear <- Nfq 
        gear_freq <- 1:Nfq
      } else {
        Ngear <- 1 
        gear_freq <- rep(1, Nfq)
        Catch <- rep(1, Nfq)
      }
    }

    if (is.null(gear_names)) {
      gear_names <- paste0("Gear_", as.character(1:Ngear))
    } else
      if ((length(gear_names) != Ngear) | !is.character(gear_names)) {
        stop("Error: gear_names must be a character vector with the same length as the number of gears.\n")
      }
    
    # Parse Catches
    if (Ngear==1) {
      catch_prop <- rep(1, Nfq)
      Ffq <- seq(Nfq)
    } else {
      if (!is.vector(Catch) | length(Catch) != Nfq | any(Catch<0) | !any(Catch>0)) {
        stop("Error: Relative catches `Catch`, including zeros where negligible, must be provided for each length frequency. \n")
      }
      sums <- tapply(Catch, period_freq, sum)
      catch_prop <- Catch[Catch > 0] / sums[period_freq[Catch > 0]] # Normalise
      names(catch_prop) <- NULL
      # Create the F index based on catches
      Ffq <- integer(Nfq)
      fi <- 1L
      for (qi in 1:Nfq) {
        if (Catch[qi] > 0) {
          Ffq[qi] <- fi
          fi <- fi + 1L
        }
      }
    }
    
    # Frequency data names
    if (is.null(freq_names)) {
      freq_names <- dplyr::case_when(
        Nperiod==1 ~ gear_names[gear_freq],
        Ngear==1   ~ period_names[period_freq],
        TRUE       ~ paste(gear_names[gear_freq], 
                           period_names[period_freq]))
    } else
      if ((length(freq_names) != Nfq) | !is.character(freq_names)) {
        stop("Error: freq_names must be a character vector with the same length as the number of data frequencies.\n")
    }

    dl <- list(
      model_name = model_name,
      # Number of gears: separate length frequencies
      NQ = Nfq,
      NG = Ngear,
      NS = 0,
      NT = Nperiod,
      # Number of F's: gears associated with non-zero catches
      NF = length(catch_prop),
      # Number of length bins
      NB = NB,
      # Place holders for number of selectivity parameters
      NP = 0,
      NM = 0,
      # Data frequency names
      fqname = as.array(freq_names),
      # Gear names
      gname = as.array(gear_names),
      # Time period names
      tpname = as.array(period_names),
      # Lower length boundaries for each bin
      LLB = LLB,
      # Length bin mid points (for plotting etc.)
      LMP = LMP,
      # Frequency data with time period and gear indexing
      fq = fq,
      Gi = as.array(gear_freq),
      Ti = as.array(period_freq),
      # Estimated total relative catch in numbers of fish, excluding zeroes for surveys etc.
      prop_catch = as.array(catch_prop),
      # Index of the Fk associated with each frequency: 0 implies catch negligible
      Fkq = as.array(Ffq),
      
      fSel = integer(0),
      GSbase = as.array(integer(Ngear)),
      GSmix1 = integer(Ngear*2),
      GSmix2 = integer(0)
    )

    dl <- blip_Linf(dl, Linf)
    dl <- blip_Galpha(dl, c(log(1 / 0.1 ^ 2), 0.25))
    dl <- blip_LH(dl, a, b, L50, L95, ma_L, wt_L, set_defaults=TRUE)
    
    if (is.null(Mk)) {
      # from Prince et al. 2015
      Mk <-
        with(dl, b * (1 - (L50 / poLinfm)) / (L50 / poLinfm))
      warning(paste0("Default Mk, based on life history invariant estimate, is: ",
                     format(Mk, digits=2), "\n"))
    } else {
      # Natural mortality
      if (Mk <= 0) 
        stop("Error: Mk must be greater than zero. \n")
    }

    dl <- blip_Mk(dl, c(log(Mk), 0.1), ref_length)
    dl <- blip_Fk(dl, NULL, 2.0)  # loose prior for fully exploited stock

    # Selectivity functions
    dl <- blicc_selfun(dl, sel_fun)
    dl <- blicc_gear_sel(dl, gear_sel)

    # Negative binomial lognormal mean
    # Medium level of overdispersion
    dl <- blip_NBphi(dl, c(log(100), 0.5))
    #  Catch: sigma for lognormal
    dl$polCs <- 0.10

    dl <- blip_sel_auto(dl)
    
    # Gauss-Laguerre quadrature grid size
    if (is.null(NK)) {
      if (any(is.na(dl$polSm)))
        dl$NK <- 110  # safe value
      else {
        df <- with(dl,
                   tidyr::expand_grid(
                     Linf = c(poLinfm - 3.1 * poLinfs, poLinfm, poLinfm + 3.1 * poLinfs),
                     Galpha = exp(c(
                       polGam - 3.1 * polGas, polGam, polGam + 3.1 * polGas
                     )),
                     Mk = exp(c(
                       polMkm - 3.1 * polMks, polMkm, polMkm + 3.1 * polMks
                     )),
                     `.draw` = 0
                   ))
        df$Gbeta <- df$Galpha / df$Linf
        df <-
          rbind(
            cbind(df, tibble::tibble(Fk = list(
              exp(dl$polFkm - 3.1 * dl$polFks)
            ))),
            cbind(df, tibble::tibble(Fk = list(exp(
              dl$polFkm
            )))),
            cbind(df, tibble::tibble(Fk = list(
              exp(dl$polFkm + 3.1 * dl$polFks)
            )))
          )
        
        df <-
          rbind(
            cbind(df, tibble::tibble(Sm = list(
              exp(dl$polSm - 3.1 * dl$polSs)
            ))),
            cbind(df, tibble::tibble(Sm = list(exp(
              dl$polSm
            )))),
            cbind(df, tibble::tibble(Sm = list(
              exp(dl$polSm + 3.1 * dl$polSs)
            )))
          )
        # Gauss-Laguerre rule
        dl$NK <- LG_Nodes(dl, df)
      }
    } else {
      dl$NK <- NK
      if (NK < 50)
        warning("Having fewer than 50 nodes for the Gauss Laguerre quadarture rule is not advised. \n")
    }
    glq <-
      statmod::gauss.quad(dl$NK, kind = "laguerre", alpha = 0.0)
    dl$gl_nodes <- glq$nodes
    dl$gl_weights <- glq$weights
    return(dl)
  }


#' Sets data object from [blicc_dat] to have new selectivity functions
#'
#' Checks selectivity functions are valid, then updates the data object with the
#' new functions so that the parameters and parameter references for the
#' selectivity functions are valid. New references are replaced in the data
#' object, which is then returned.
#'
#' Particular selectivity functions can be changed using the `sel_indx` parameter,
#' which indexes the relevant functions. In this case, only those
#' functions will need to have the prior parameters updated.
#'
#' @export
#' @inheritParams blicc_mpd
#' @inheritParams blicc_dat
#' @param sel_indx  An integer vector index of the selectivity functions to be
#'   changed. If not provided, all functions will be replaced.
#' @return The data object `blicc_ld` with the new selectivities
#' @examples
#' new_ld <- blicc_selfun(gillnet_ld, sel_fun="logistic", model_name = "Logistic Selectivity")
#'
blicc_selfun <-
  function(blicc_ld,
           sel_fun,
           sel_indx = NULL,
           model_name = NULL) {

    sel_fun <- parse_selectivity(sel_fun, blicc_ld)
    mixwt <- NULL
    mixwts <- NULL
    
    if (is.null(sel_indx)) {
      nfSel <- sel_fun
      if (blicc_ld$NM > 0) {
        warning("Gear selectivity mixtures have been removed \n")
        blicc_ld$NM <- 0
        GSmix1 <- integer(blicc_ld$NG * 2)
        GSmix2 <- integer(0)
      }
    } else {
      if (blicc_ld$NS == 0)
        stop(
          "Error: No previous functions to reference, so replace all functions (sel_indx=NULL). \n"
        )
      sel_indx <- parse_sel_indx(sel_indx, blicc_ld, FALSE)
      if (length(sel_indx) != length(sel_fun))
        stop("Error: sel_indx must have the same length as sel_fun. \n")
      
      nfSel <- blicc_ld$fSel
      nfSel[sel_indx] <- sel_fun
      if (blicc_ld$NM > 0) {
        # mixture weights
        mixwt <-
          blicc_ld$polSm[(blicc_ld$NP + 1L):length(blicc_ld$polSm)]
        mixwts <-
          blicc_ld$polSs[(blicc_ld$NP + 1L):length(blicc_ld$polSm)]
      }
    }
    # copy parameters
    NS <- length(nfSel)
    npar <- Rsel_functions()$npar[nfSel]
    NP <- sum(npar)
    nSm <- rep(NA_real_, NP) # lognormal mu parameters
    nSs <- nSm
    spari <- integer(NS)
    spare <- spari
    np <- 1L
    for (i in 1:NS) {
      spari[i] <- np
      np <- np + npar[i]
      spare[i] <- np - 1L
      if (i <= blicc_ld$NS) {
        if (nfSel[i] == blicc_ld$fSel[i]) {
          # copy parameters
          nSm[spari[i]:spare[i]] <-
            with(blicc_ld, polSm[sp_i[i]:sp_e[i]])
          nSs[spari[i]:spare[i]] <-
            with(blicc_ld, polSs[sp_i[i]:sp_e[i]])
        }
      }
    }
    # Replace function list
    blicc_ld$fSel <- as.array(nfSel)
    blicc_ld$NS <- NS
    blicc_ld$sp_i <- as.array(spari)    #start
    blicc_ld$sp_e <- as.array(spare)    #end
    blicc_ld$NP <- NP
    blicc_ld$polSm <- c(nSm, mixwt)
    blicc_ld$polSs <- c(nSs, mixwts)
    if (!is.null(model_name))
      blicc_ld$model_name <- model_name
    
    if (all(blicc_ld$GSbase > 0))
      blicc_ld <- blip_sel_auto(blicc_ld, sel_indx = sel_indx)
    return(blicc_ld)
  }


#' Returns data list from [blicc_dat] with new gear-selectivity links
#'
#' Sets up the data list with the new links between the selectivity functions
#' and the gears. This allows the gears to be associated with more than one
#' selectivity function (mixtures), so increasing the flexibility of the gear
#' selectivity. The new references are replaced in the data object, which is
#' then returned.
#'
#' @export
#' @inheritParams blicc_dat
#' @inheritParams blicc_mpd
#' @param gear  The names or integer vector of the gears being linked to the
#'   selectivity functions
#' @return The data object blicc_ld with the new selectivity-gear links.
#'
blicc_gear_sel <-
  function(blicc_ld,
           gear_sel = NULL,
           gear = NULL,
           model_name = NULL) {

    ErrorMsg <- paste0("Error: the specified selectivity function must equal the number of gears (",
                       length(gear), "or a list of integer vectors must be provided to gear_sel to link gears to selectivity. \n")

    if (blicc_ld$NS==0) stop("Error: Selectivity functions must be set first. \n")

    if (is.null(gear_sel)) {
      if (blicc_ld$NS != blicc_ld$NG)
        stop("Error: if gear_sel is not provided, the number of selectivity functions must equal the number of gears. \n")

      # 1 selectivity function for each gear
      if (! is.null(gear)) warning("gear parameter is ignored. \n")
      blicc_ld$GSbase <- as.array(1L:blicc_ld$NS)
      gsmix1 <- integer(2L*blicc_ld$NG)
      gsmix2 <- integer(0)
    } else {
      if (! (is.list(gear_sel) & is.vector(gear_sel)))
        stop("Error: gear_sel must be a list of integer vectors linking the gear to the relevant selectivity functions. \n")
      if (is.null(gear)) {
        if (length(gear_sel) != blicc_ld$NG)
          stop("Error: If gear parameter is not given, gear_sel must be the same length as the number of gears. \n")
        gear <- 1L:blicc_ld$NG
      } else {
        gear <- parse_gear(gear, blicc_ld)
        if (length(gear_sel) != length(gear))
          stop("Error: the length of the gear parameter must equal the length of gear_sel. \n")
      }
      # Mixtures may be present

      if (is.null(blicc_ld$GSbase)) {
        # Initialise gear links
        if (length(gear) != blicc_ld$NG)
          stop("Error: Gear links must be initialised for all gears. \n")
        blicc_ld$GSbase <- array(0L, blicc_ld$NG)
        blicc_ld$GSmix1 <- integer(2L*blicc_ld$NG)
        blicc_ld$GSmix2 <- integer(0)
        blicc_ld$NM <- 0
      }

      # Rebuild links
      gsmix1 <- integer(2L*blicc_ld$NG)
      gsmix2 <- integer(0)
      for (gi in 1L:blicc_ld$NG) {
        # remove any mixtures for that gear
        si <- 1L + 2L*(gi-1L)
        if (gi %in% gear) {
          #replace
          gii <- which(gear==gi)
          new_map <- gear_sel[[gii]]
          if (! (is.vector(new_map) & is.numeric(new_map)))
            stop("Error: New links to selectivity must be provided as integer vectors indexing selectivity functions. \n")
          new_map <- unique(round(new_map))
          if ((length(new_map) < 1L) | (min(new_map) < 1L) | (max(new_map) > blicc_ld$NS))
            stop(paste0("Error: Selectivity function index out of bounds for gear ", as.character(gii), " \n"))
          blicc_ld$GSbase[gi] <- new_map[1L]
          if (length(new_map) > 1L) {
            # mixtures present
            n <- length(gsmix2) + 1L
            L <- length(new_map) - 2L
            gsmix2 <- c(gsmix2, new_map[-1])
            gsmix1[si] <- n
            gsmix1[si+1] <- n + L
          }
        } else {
          # keep current links
          if (blicc_ld$GSmix1[si] > 0) {
            n <- length(gsmix2) + 1L
            L <- blicc_ld$GSmix1[si+1] - blicc_ld$GSmix1[si]
            gsmix2 <- c(gsmix2, blicc_ld$GSmix2[blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]])
            gsmix1[si] <- n
            gsmix1[si+1] <- n + L
          }
        }
      }
    }
    # replace indexes
    blicc_ld$GSmix1 <- gsmix1
    blicc_ld$GSmix2 <- as.array(gsmix2)
    blicc_ld$NM <- length(gsmix2)
    # extend parameters for mixtures
    blicc_ld$polSm <- c(blicc_ld$polSm[1:blicc_ld$NP], rep(NA_real_, blicc_ld$NM))
    blicc_ld$polSs <-  c(blicc_ld$polSs[1:blicc_ld$NP], rep(1.5, blicc_ld$NM))

    if (!is.null(model_name))
      blicc_ld$model_name <- model_name

    if ( blicc_model_OK(blicc_ld) == "OK" )
      blicc_ld <- blip_sel_auto(blicc_ld)

    return(blicc_ld)
  }



#' Returns data list from [blicc_dat] with only the selected time period
#' included
#'
#' Subsets the data list with only selectivity components, gears, frequencies
#' and parameters relevant to the defined time period.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param time_period  The names or integer vector of the time periods being
#'   selected
#' @return The data object blicc_ld subset for the new time period.
#' 
blicc_period_filter <- function(blicc_ld, time_period) {
  time_period <- parse_period(sort(time_period), blicc_ld)
  
  if (blicc_ld$NT==1L) return(blicc_ld)
  
  Qindx <- blicc_ld$Ti %in% time_period
  Findx <- Qindx & blicc_ld$Fkq > 0
  Gindx <- unique(blicc_ld$Gi[Qindx])
  Sindx <- get_selectivities(Gindx, blicc_ld)
  Pindx <- integer(0)
  for (si in Sindx) {
    Pindx <- c(Pindx, blicc_ld$sp_i[si]:blicc_ld$sp_e[si])
  }
  #Mindx <- + MIXWeights
  
  ld <- blicc_ld
  ld$model_name <- paste(blicc_ld$model_name, "Subset Periods", 
                         paste(as.character(time_period), collapse=" "))
  ld$NQ <- sum(Qindx)
  ld$NG <- length(Gindx)
  ld$NS <- length(Sindx)
  ld$fSel <- as.array(blicc_ld$fSel[Sindx])

  ld$NT <- length(time_period)
  ld$NF <- sum(Findx)
  ld$NP <- length(Pindx)
  #ld$NM <- length(Mindx)
  ld$fqname <- as.array(blicc_ld$fqname[Qindx])
  ld$gname <- as.array(blicc_ld$gname[Gindx])
  ld$tpname <- as.array(blicc_ld$tpname[time_period])
  
  ld$fq <- blicc_ld$fq[Qindx]
  ld$Gi <- as.array(match(blicc_ld$Gi[Qindx], Gindx))
  ld$Ti <- as.array(match(blicc_ld$Ti[Qindx], time_period))
  
  ld$prop_catch <- as.array(blicc_ld$prop_catch[Qindx])
  ld$Fkq <- as.array(match(blicc_ld$Fkq[Qindx], which(Findx)))
  ld$Fkq[is.na(ld$Fkq)] <- 0
  sums <- with(ld, tapply(prop_catch, Ti[Fkq>0], sum))
  ld$prop_catch <- with(ld, prop_catch / sums[Ti[Fkq>0]]) # Normalise
  
  
  #seq_along(Gindx)
  ld$GSbase <- as.array(blicc_ld$GSbase[Gindx]) # resequence
  mxn <- integer(0)
  mxpar <- mxpars <- double(0)
  
  ld$GSmix1 <- integer(2*length(Gindx))
  mix_1 <- 1L
  mix_2 <- 1L
  for (gi in seq(Gindx)) {
    mgi <- (gi-1L)*2L + 1L
    if (blicc_ld$GSmix1[mgi] > 0) {
      mxindx <- blicc_ld$GSmix2[blicc_ld$GSmix1[gi]:blicc_ld$GSmix1[gi+1L]]
      ld$GSmix1[mix_1] <- mix_2
      mix_1 <- mix_1 + 1L
      mix_2 <- mix_2 + length(mxindx) - 1L
      ld$GSmix1[mix_1] <- mix_2
      mix_1 <- mix_1 + 1L
      mix_2 <- mix_2 + 1L
      mxn <- c(mxn, match(mxindx, Sindx))
      mxpar <- c(mxpar, blicc_ld$polSm[blicc_ld$NP+blicc_ld$GSmix1[mgi]:blicc_ld$GSmix1[mgi+1L]])
      mxpars <- c(mxpars, blicc_ld$polSs[blicc_ld$NP+blicc_ld$GSmix1[mgi]:blicc_ld$GSmix1[mgi+1L]])
    } else {
      mix_1 <- mix_1 + 2L
    }
    
  }
  ld$GSmix2 <- mxn
  ld$NM <- length(mxn)
  ld$polFkm <- as.array(blicc_ld$polFkm[Findx])

  npar <- Rsel_functions()$npar[ld$fSel]
  spare <- spari <- integer(ld$NS)
  np <- 1L
  for (i in seq_along(ld$fSel)) {
    spari[i] <- np
    np <- np + npar[i]
    spare[i] <- np - 1L
  }
  # Replace function list
  ld$sp_i <- as.array(spari)    #start
  ld$sp_e <- as.array(spare)    #end

  ld$polSm <- c(blicc_ld$polSm[Pindx], mxpar)
  ld$polSs <- c(blicc_ld$polSs[Pindx], mxpars)
  
  return(ld)    
}

#' Default start parameters centred on priors
#'
#' Provides a list of start parameters for the MCMC - defaults to zero values
#'
#' @inheritParams blicc_mpd
#' @return A list of initial parameter values for the BLICC Stan model
#' @noRd
#'
blicc_ini <- function(blicc_ld) {
  return(list(
    nLinf   = 0,
    nGalpha = 0,
    nMk    = 0,
    nFk    = as.array(rep(0.0, blicc_ld$NF), dim=1),
    nSm    = rep(0.0, blicc_ld$NP+blicc_ld$NM),
    nNB_phi = 0
  ))
}


#' Provides start parameters to the number of chains required for the MCMC run
#'
#' Provides a list of start values for the BLICC model MCMC using function
#' `blicc_ini()` if no parameters are provided.
#'
#' @inheritParams blicc_mpd
#' @param  pchain A number indicating the number of chains. Default is 4
#' @param  par    A list of start parameter values (optional)
#' @return A list of start parameters equal to the number of chains for the
#'   BLICC Stan model
#' @noRd
#'
blicc_mcmc_ini <- function(blicc_ld, pchain = 4, par = NULL) {
  if (is.null(par)) {
    return(rep(list(blicc_ini(blicc_ld)), pchain))
  } else {
    return(rep(list(par), pchain))
  }
}


#' Conducts basic checks of the data object from [blicc_dat]
#'
#' The function returns "OK" or an error message. Checks are not exhaustive.
#'
#' @inheritParams blicc_mpd
#' @return Text reporting whether data object is OK, and if not why.
#' @noRd
#'
blicc_model_OK <- function(blicc_ld) {
  if (is.null(blicc_ld$gl_nodes))
    return("Error: Gauss Laguerre nodes missing. \n")
  rsel <- Rsel_functions()
  # Check indexing
  if (any(blicc_ld$fSel > length(rsel$npar)) | any(blicc_ld$fSel <= 0))
    return("Error: Selectivity function(s) invalid. \n")

  fref <- c(blicc_ld$GSbase, blicc_ld$GSmix2)
  if (any(fref > blicc_ld$NS) | any(fref <= 0))
    return("Error: Gear selectivity reference is out of range. \n")
  if (any(is.na(match(1:blicc_ld$NS, fref))))
    return("Error: A selectivity function is not referenced. \n")
  if (any(is.na(blicc_ld$polSm)))
    return("Warning: Selectivity priors are not set. \n")

  return("OK")
}



