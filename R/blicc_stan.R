# Model Fit Functions -----------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

#' Finds Bayesian length interval catch curve maximum posterior density point
#'
#' The maximum posterior density point is found using on the likelihood for the
#' length frequency data and the priors.  The model assumes constant
#' recruitment, double-sided normal selectivity, a single natural mortality
#' parameter, and Gamma-distributed von Bertalanffy growth to define expected
#' numbers of fish in pre-defined length bins. The model also assumes
#' a negative binomial likelihood function to fit to a single
#' length frequency sample. Selectivity, mortality, asymptotic mean length and
#' scale parameters are fitted to generate the spawning potential ratio as a
#' measure of the stock status.
#'
#' It is recommended to test the model and data object using this procedure
#' before investing in the MCMC fit using [blicc_fit()].
#'
#' @details
#' Unlike [blicc_fit()], the model finds a point estimate. This is much faster
#' than the MCMC, but does not provide full information on uncertainty.
#' The model estimates mortality and survival through sequential length
#' intervals. This is then used to derive abundance and catch in numbers within
#' each length bin. The mortality is required to remain constant within each
#' bin, but otherwise can vary arbitrarily. However, in this implementation
#' mortality is constrained by a parametric model with constant natural
#' mortality by length and a double-sided normal selectivity function for
#' fishing mortality. See [blicc_dat()] and [blicc_fit()] documentation for
#' a description of the model and data requirements.
#'
#' The standard errors are estimated from 1000 random draws from a multivariate
#' normal with mean at the MPD mode and covariance estimated from the inverted
#' Hessian matrix. This is an approximation and SE estimates may differ from
#' the full MCMC sampling.
#'
#' In addtion to the parameters, the fit reports the spawning potential ratio
#' (SPR) and the log probability at the mode (`lp__`).
#'
#' @export
#' @param  blicc_ld    A standard data list suitable for the model (see function
#'   [blicc_dat])
#' @return A tibble of parameter estimates (mpd) with standard errors.
#' @examples
#' mpd_fit <- blicc_mpd(eg_ld)
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
                      paste0("Sm[", as.character(1:blicc_ld$NP), "]"),
                      "NB_phi", "Gbeta", "SPR")
  rd <- fit$theta_tilde[ , full_par_names]
  if (blicc_ld$NF>1)
    par_names <- c("Linf", "Galpha", "Mk",
                   paste0("Fk", as.character(1:blicc_ld$NF)),
                   paste0("Sm", as.character(1:blicc_ld$NP)),
                   "NB_phi", "Gbeta", "SPR")
  else
    par_names <- c("Linf", "Galpha", "Mk", "Fk",
                   paste0("Sm", as.character(1:blicc_ld$NP)),
                   "NB_phi", "Gbeta", "SPR")
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
#' The model is fitted to length frequency data using a Markov chain Monte
#' Carlo (MCMC) simulation. The model assumes constant recruitment,
#' double-sided normal selectivity, a single natural mortality parameter, and
#' Gamma-distributed von Bertalanffy growth to define expected numbers of fish
#' in pre-defined length bins. The model also assumes a negative binomial
#' likelihood function to fit to a single length frequency sample. Selectivity,
#' mortality, asymptotic mean length and scale parameters are fitted to
#' generate the spawning potential ratio as a measure of the stock status.
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
#' The model fits 4 parameters plus 3 or more parameters for each selectivity
#' function being fitted. In general, there is insufficient support from
#' a single length frequency samples to estimate all these with any precision,
#' and therefore informative priors are used where necessary. The default
#' values for priors can follow recommendations based on meta-analyses or
#' life-history invariant (see [blicc_dat]). However the user is still
#' required to provide a prior on the mean maximum length (Linf) in each case.
#'
#' It is recommended to fit a model using the [blicc_mpd()] function before
#' attempting the MCMC fit. This should show up any problems with the model
#' or data object.
#'
#' The fit starts at the posterior mode which is found using the Stan optimizer.
#' Then a default 4 MCMC chains are initiated and run concurrently. Comparing
#' the chains provides the diagnostics on whether the MCMC has converged or not.
#'
#' It is quite possible to get messages of the type
#' "Exception: neg_binomial_2_lpmf: ..." when the fitting starts. This should
#' not be a problem as long as the message does not persist. If it does either
#' there is a problem with the input data or the MCMC settings - particularly
#' the priors. The priors are setup in the [blicc_dat] function. The MCMC
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
#'
#' For example, the parameter "control=list(adapt_delta=0.9)" can be added
#' to the parameter list on calling the function.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param  ntarget     Target draw for the MCMC
#' @param  nwarmup     Warm up iterations for the Stan MCMC
#' @param  nchain      Number of chains to run in parallel
#' @param  ...         Other arguments passed to `rstan::sampling()`.
#' @return An object of class stanfit returned by `rstan::sampling()`
#' @examples
#' stf <- blicc_fit(eg_ld, ntarget=100, nwarmup=200, nchain=1) #Quick test
#'
blicc_fit <- function(blicc_ld,
                      ntarget = 2000,
                      nwarmup = 1000,
                      nchain = 4,
                      ...) {

  options(mc.cores = parallel::detectCores())  # Needed for parallel chains

  #### remove for the package version
  # Test version
  # stan_fn <- here("R", "BLICC_mg.stan")
  # stmod <- stan_model(stan_fn, model_name = "BLICC")
  #

  # Find the posterior mode to start
  # # Distributed version
  res <-
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
    init = blicc_mcmc_ini(blicc_ld, nchain, res$par),
    verbose = FALSE,
    ...
  )
  return(stf)
}


#' Generates a data list and initial parameter values for the BLICC model
#'
#' Vectors for the lower bound for each length bin and the number of fish in
#' each length bin and combined with prior parameters into a data list as
#' expected by the BLICC model.
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
#' The length frequencies should include zeros for bins which
#' contained no fish. The length frequency should be bounded by a single zero
#' bin as the first and last bin in the frequency.
#'
#' To allow maximum flexibility, prior hyper-parameters for the selectivity
#' functions (as log-values for means in the log-normal distribution) are
#' defined in a single vector. These are linked to each function using an
#' integer matrix with a row for each selectivity function (gear)
#'
#' @export
#' @param model_name A name for the model (species or fishery). Optional.
#' @param LLB  A vector of lower length boundaries for each length frequency
#'   bin. Required.
#' @param fq   A list of vectors, each the same length as LLB, containing the
#'   frequency data. Zeroes must be included. If only one gear, can be a vector.
#'   Required.
#' @param Linf A vector of two values: mean and sd of the prior maximum mean
#'   length for the stock. Required.
#' @param sel_fun A vector of selectivity function names (character) or index
#'   (1-4), with 1="logistic", 2="normal", 3="ss_normal", 4="ds_normal" for each
#'   vector in fq. Required.
#' @param Catch A vector of relative catches. Required, if the number of gears
#'   is more than one.
#' @param gear_names A vector of names for the gears for reference. Optional.
#' @param Mk   Natural mortality divided by the growth rate K (usually around
#'   1.5). Optional.
#' @param ref_length  The reference length for the inverse length function if it
#'   is used. (default is NA i.e. fixed natural mortality)
#' @param wt_L  Weight (biomass) for each length bin. Optional.
#' @param ma_L  Mature biomass for each length bin. Optional.
#' @param a  The length-weight parameter: a*L^b. Not used in model fitting, but
#'   for plots etc. Optional.
#' @param b  The length-weight exponent (a*L^b, usually close to 3.0) Optional.
#' @param L50 The length at 50% maturity, often referred to as "first" maturity,
#'   and primarily applies to females. Must be less than Linf, usually around
#'   0.66 Linf, which is assumed if it is not provided. Optional.
#' @param L95  The length at 95% maturity. Must be greater than L50 and less
#'   than Linf. A small increment is added L50 if it is not provided. Optional.
#' @param NK   Number of nodes for the Gauss Laguerre quadrature rule. 110 is a
#'   safe value, but extends the run time of all calculations.  Optional.
#' @return     A list of vectors and numbers structured suitable for use in the
#'   Stan model `BLICC.stan`.
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
           Catch = 1,
           gear_names = NA,
           Mk = NA,
           ref_length = -1,
           wt_L = NA,
           ma_L = NA,
           a = 1,
           b = 3,
           L50 = NA,
           L95 = NA,
           NK = NA) {

    if (!is.vector(LLB, mode="numeric")) {
      stop(
        "Error: please supply the lower bound of bins as a vector of unique values and in ascending order."
      )
    }
    if (any(diff(LLB) <= 0)) {
      stop("Error: lower bound of bins must be unique and in ascending order.")
    }
    NB <- length(LLB)
    LMP <- c((LLB[-NB] + LLB[-1]) * 0.5, LLB[NB] + 0.5*(LLB[NB]-LLB[NB-1])) # Length bin mid points (LMP) used for plotting etc. Not used in Stan model.

        # Data vectors
    if (is.vector(fq, mode="numeric")) fq <- list(fq)  #convert to list if possible
    if (!is.list(fq)) {
      stop(
        "Error: frequency data must be supplied as a list of vectors (or can be given as a vector if only one.)"
      )
    }
    Ngear <- length(fq)
    for (i in 1:Ngear)
      if (!is.vector(fq[[i]], mode="numeric") | (NB != length(fq[[i]]))) {
        stop(
          "Error: please provide each frequency in each bin, including zeroes as vectors of the same length as the number of bins."
        )
      }

    # Catches
    if (!is.vector(Catch) | length(Catch) != Ngear | any(Catch<0) | !any(Catch>0)) {
      stop("Error: Relative catches must be provided for each length frequency.")
    }
    Fgear <- integer(Ngear)
    fi <- 1
    for (gi in 1:Ngear) {
      if (Catch[gi] > 0) {
        Fgear[gi] <- fi
        fi <- fi + 1
      }
    }

    catch_prop <- Catch[Catch > 0]/sum(Catch)

    # Gear names
    if (any(is.na(gear_names))) {
      gear_names <- paste0("Gear_", as.character(1:Ngear))
    } else
      if ((length(gear_names) != Ngear) | !is.character(gear_names)) {
        stop("Error: gear_names must be a character vector with the same length as the number of gears.")
      }

    # Natural mortality
    if (!is.na(Mk) & Mk <= 0) {
      stop("Error: Mk must be greater than zero.")
    }

    dl <- list(
      model_name = model_name,
      # Number of gears: separate length frequencies
      NG = Ngear,
      # Number of F's: gears associated with non-zero catches
      NF = length(catch_prop),
      # Number of length bins
      NB = NB,
      # Gear names
      gname = as.array(gear_names),
      # Lower length boundaries for each bin
      LLB = LLB,
      # Length bin mid points (for plotting etc.)
      LMP = LMP,
      # Frequency data
      fq = fq,
      # Estimated total relative catch in numbers of fish, excluding zeroes for surveys etc.
      prop_catch = as.array(catch_prop),
      # Index of the Fk associated with each gear: 0 implies catch negligible
      Fkg = as.array(Fgear)
    )
    # Selectivity functions
    dl <- blicc_selfun(dl, 1:Ngear, sel_fun)
    dl <- blip_Linf(dl, Linf)
    dl <- blip_LH_param(dl, a, b, L50, L95, ma_L, wt_L, set_defaults=TRUE)
    dl <- blip_Galpha(dl, c(log(1 / 0.1 ^ 2), 0.25))

    if (is.na(Mk)) {
      # from Prince et al. 2015
      Mk <-
        with(dl, b * (1 - (L50 / poLinfm)) / (L50 / poLinfm))
      warning(paste0("Warning: Default Mk, based on life history invariant estimate, is: ",
                     format(Mk, digits=2)))
    }

    dl <- blip_Mk(dl, c(log(Mk), 0.1), ref_length)
    dl <- blip_Fk(dl, NA, 2.0)  # loose prior for fully exploited stock
    dl <- blip_selectivity(dl)

    # Negative binomial lognormal mean
    # Medium level of overdispersion
    dl <- blip_NBphi(dl, c(log(100), 0.5))
    #  Catch: sigma for lognormal
    dl$polCs <- 0.01

    # Gauss-Laguerre quadrature grid size
    if (is.na(NK)) {
      df <- with(dl,
                 tidyr::expand_grid(
                   Linf = c(poLinfm-3.1*poLinfs, poLinfm, poLinfm+3.1*poLinfs),
                   Galpha = exp(c(polGam-3.1*polGas, polGam, polGam+3.1*polGas)),
                   Mk = exp(c(polMkm-3.1*polMks, polMkm, polMkm+3.1*polMks)),
                   `.draw` = 0
                 ))
      df$Gbeta <- df$Galpha/df$Linf
      df <- rbind(cbind(df, tibble::tibble(Fk=list(exp(dl$polFkm-3.1*dl$polFks)))),
                  cbind(df, tibble::tibble(Fk=list(exp(dl$polFkm)))),
                  cbind(df, tibble::tibble(Fk=list(exp(dl$polFkm+3.1*dl$polFks)))))

      df <- rbind(cbind(df, tibble::tibble(Sm=list(exp(dl$polSm-3.1*dl$polSs)))),
                  cbind(df, tibble::tibble(Sm=list(exp(dl$polSm)))),
                  cbind(df, tibble::tibble(Sm=list(exp(dl$polSm+3.1*dl$polSs)))))
      # Gauss-Laguerre rule
      dl$NK <- LG_Nodes(dl, df)
    } else {
      dl$NK <- NK
      if (NK < 50)
        warning(
          "Warning: Having fewer than 50 nodes for the Gauss Laguerre quadarture rule is not advised."
        )
    }

    glq <-
      statmod::gauss.quad(dl$NK, kind = "laguerre", alpha = 0.0)
    dl$gl_nodes <- glq$nodes
    dl$gl_weights <- glq$weights
    return(dl)
  }


#' Sets data object from [blicc_dat()] to have a new selectivity functions
#'
#' Checks the gear and selectivity are valid, then sets up the data object
#' so that the parameters and parameter references for the selectivity
#' functions are valid. New references are replaced in the data object,
#' which is then returned.
#'
#' @export
#' @inheritParams blip_Linf
#' @param Gear     The gears for which the selectivity functions are being
#'   changed
#' @param sel_fun  The selectivity functions that gears are being set to either
#'   as integers or strings.
#' @return The data object blicc_ld with the new selectivities
#' @examples
#' new_ld <- blicc_selfun(eg_ld, Gear=1, sel_fun="logistic", model_name = "Logistic Selectivity")
#'
blicc_selfun <-
  function(blicc_ld,
           Gear,
           sel_fun,
           model_name = NA) {
    Gear <- parse_gear(Gear, blicc_ld)
    sel_fun <- parse_selectivity(sel_fun, blicc_ld)
    if (length(sel_fun) != length(Gear)) {
      stop(paste0("Error: the specified selectivity function must equal the number of gears (",
                  length(Gear)))
    }
    if (is.null(blicc_ld$fSel)) blicc_ld$fSel <- as.array(integer(blicc_ld$NG))

    blicc_ld$fSel[Gear] <- sel_fun

    npar <- Rsel_functions()$npar[blicc_ld$fSel]

    blicc_ld$NP <- sum(npar)
    spar <- integer(blicc_ld$NG)
    spare <- spar
    np <- 1L
    for (i in 1:blicc_ld$NG) {
      spar[i] <- np
      np <- np + npar[i]
      spare[i] <- np-1L
    }
    blicc_ld$sp_i <- spar    #start
    blicc_ld$sp_e <- spare   #end
    if (!is.na(model_name))
      blicc_ld$model_name <- model_name
    blicc_ld <- blip_selectivity(blicc_ld)
    return(blicc_ld)
  }


#' Sets data object from [blicc_dat()] to have a new `Linf` prior
#'
#' Checks the new `Linf` is a vector of 2 values, the new mean and sd for
#' the normal prior to be used. These are replaced in the data object,
#' which is then returned.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param Linf  A numeric vector of double containing the mean and
#' sd for the prior normal.
#' @param model_name A string for a replacement model name in the data
#' object. Optional.
#' @return The data object blicc_ld but with the prior for `Linf` changed.
#' @examples
#' new_ld <- blip_Linf(eg_ld, c(30,2), model_name="Sensitivity")
#'
blip_Linf <- function(blicc_ld,
                      Linf,
                      model_name=NA) {
  if (! is.vector(Linf, mode="double") | length(Linf) != 2)
    stop("Error: Linf must be a vector of 2 (mu, sigma) for the prior.")
  if (Linf[1] <= min(blicc_ld$LLB)) {
    stop(paste(
      "Error: Linf must be greater than the lowest bin value:",
      as.character(min(blicc_ld$LLB))
    ))
  }
  if (Linf[2] <= 0) {
    stop("Error: Linf prior sd must be greater than zero.")
  }

  if (!is.na(model_name))
    blicc_ld$model_name <- model_name
  # Expected Linf
  blicc_ld$poLinfm <- Linf[1]
  # sd for the normal Linf, see above
  blicc_ld$poLinfs <- Linf[2]
  return(blicc_ld)
}


#' Sets data object from [blicc_dat()] to have a new Galpha prior
#'
#' The Galpha prior is updated with new values.
#' CV: could be 5% (log(1/0.05^2)) to 30% (log(1/0.3^2)).
#' A 30% CV makes length very uninformative on age however.
#' Default CV=10%
#'
#' @export
#' @inheritParams blip_Linf
#' @param lGalpha  A numeric vector of double containing the mean and sd for the
#'   prior log-normal.
#' @return The data object blicc_ld but with the prior for Galpha changed.
#'
blip_Galpha <- function(blicc_ld,
                        lGalpha) {
  if (!(is.numeric(lGalpha) & length(lGalpha)==2))
    stop("Error: LGalpha must be numeric vector of the mu,sd for the Galpha lognormal prior.")
  if (lGalpha[1] > log(1/0.05^2) | lGalpha[1] < log(1/0.3^2))
    warning("The Galpha prior lognormal mean is outside the expected 5-30% CV.")
  # Growth mean CV
  blicc_ld$polGam <- lGalpha[1]
  blicc_ld$polGas <- lGalpha[2]
  return(blicc_ld)
}


#' Sets data object from [blicc_dat()] to have a new Mk prior
#'
#' Checks the new values are valid, then provide the mean and sd for
#' the log-normal prior to be used. Also, if the ref_length parameter
#' is defined, apply the inverse length function for natural mortality.
#' These are replaced in the data object, which is then returned.
#' Note that the Mk must be provided as the log value.
#'
#' @export
#' @inheritParams blip_Linf
#' @inheritParams blicc_dat
#' @param lMk   The (log) mean and sigma for the lognormal natural mortality prior
#' @return The data object blicc_ld but with the prior and function for Mk
#'   changed.
#' @examples
#' new_ld <- blip_Mk(eg_ld, lMk=c(log(1.9), 0.2), ref_length=25)
#'
blip_Mk <- function(blicc_ld,
                   lMk = as.numeric(c(NA, NA)),
                   ref_length = NA,
                   model_name = NA) {
  # Natural mortality
  if (!is.na(ref_length)) {
    if (ref_length>0){
      if (ref_length < min(blicc_ld$LLB) | ref_length > max(blicc_ld$LLB))
        stop("Error: The reference length for the natural mortality must be within the length frequencies.")
      M_L <- ref_length/blicc_ld$LMP
    } else {
      M_L <- rep(1, blicc_ld$NB)  # Fixed natural mortality
      ref_length <- -1
    }
    blicc_ld$M_L <- M_L
    blicc_ld$ref_length <- ref_length
  }

  if (! (is.vector(lMk, mode = "numeric") & length(lMk)==2))
    stop("Error: natural mortality must be provided as vector of mean and sigma for the lognormal.")

  if (!is.na(model_name))
    blicc_ld$model_name <- model_name
  if (!is.na(lMk[1])) {
    if (lMk[1]<0 | lMk[1]>log(5))
      warning(paste0("Warning: Mk outside range 1-5: ",
                     format(exp(lMk[1]), digits=2), " (make sure you provide the log-Mk)"))
    blicc_ld$polMkm <- lMk[1]
  }
  if (!is.na(lMk[2]))
    blicc_ld$polMks <- lMk[2]
  return(blicc_ld)
}


#' Sets data object from [blicc_dat()] to have a new Fk prior
#'
#' The Fk lognormal prior is updated with new mu and a single sigma parameter.
#' If no value (`NA`) is provided for lFk, the defaults are assigned. These are
#' the current natural mortality mu for fishing mortality mu and 2.0 for the
#' sigma parameter.
#'
#' @export
#' @inheritParams blip_Linf
#' @param lFk   A vector of double containing the lognormal mean Fk for each gear
#' @param lFks  A double containing the lognormal sd Fk for all gears
#' @return The data object `blicc_ld` with the prior for `Fk` changed.
#' @examples
#' new_ld <- blip_Fk(eg_ld, lFk=log(1.9), lFks=1.5)
#'
blip_Fk <- function(blicc_ld,
                   lFk  = NA,
                   lFks = NA,
                   model_name=NA) {
  if (any(is.na(lFk))) {
    lFk <- blicc_ld$polMkm + log(blicc_ld$prop_catch)
  }
  if (is.na(lFks))
    lFks = 2.0
  if (! (is.numeric(lFk) & is.numeric(lFks)))
    stop("Error in blip_Fk: supplied values are not numeric.")
  if (blicc_ld$NF != length(lFk))
    stop(paste("Error in blip_Fk: lFk array not of length: ", as.character(blicc_ld$NF)))
  if (! is.na(model_name))
    blicc_ld$model_name <- model_name
  blicc_ld$polFkm <- as.array(lFk)
  blicc_ld$polFks = lFks
  return(blicc_ld)
}


#'Sets data object from [blicc_dat()] to have a new selectivity functions
#'
#'Checks the gear and selectivity are valid, then sets up the data object so
#'that the parameters and parameter references for the selectivity functions are
#'valid. New references are replaced in the data object, which is then returned.
#'
#'@export
#'@inheritParams blip_Linf
#'@inheritParams blicc_dat
#'@param set_defaults Logical indicating whether to set defaults or not if
#'  parameters are NA
#'@return The data object blicc_ld with the new life-history parameters
#'@examples
#'new_ld <- blip_LH_param(eg_ld, a=1.2e-5, b=2.95, model_name="Alternative LW")
#'
blip_LH_param <-
  function(blicc_ld,
           a = NA,
           b = NA,
           L50 = NA,
           L95 = NA,
           ma_L = NA,
           wt_L = NA,
           model_name = NA,
           set_defaults = FALSE) {

    Linf <- blicc_ld$poLinfm
    # Weight and maturity
    if (is.na(L50)) {
      if (set_defaults)
        L50 <- 0.66 * Linf
    } else {
      if (L50 <= 0.2 * Linf | L50 >= Linf) {
        stop("Error: Length at 50% maturity must be greater than 0.2*Linf and less than Linf.")
      }
    }
    if (is.na(L95)) {
      # Ls = -log(1/0.95 - 1)/(L95-L50)
      if (set_defaults)
        Ls <- -log(1 / 0.95 - 1) / (0.05 * (Linf - L50))
      else
        Ls <- NA
    } else {
      if (L95 <= L50 | L95 >= Linf) {
        stop("Error: Length at 95% maturity must be greater than L50 and less than Linf.")
      }
      Ls <- -log(1 / 0.95 - 1) / (L95 - L50)
    }

    if (any(is.na(wt_L))) {
      if (is.na(a) & set_defaults) {
        warning("Warning: No weight-at-length information provided - the weight units will be incorrect.")
        a <- 1.0
      }
      if (is.na(b) & set_defaults)
        b <- 3.0
      else if (b <= 2 | b > 4) {
        stop("Error: Length-weight exponent (b) must be greater than 2 and less than 4.")
      }
      if (! (is.na(a) | is.na(b)))
        wt_L <- with(blicc_ld, a * exp(b*log(LMP)))    # Estimated biomass per recruit
      else
        warning("a or b not specified: weight-at-length not changed.")
    } else {
      if (length(wt_L) != blicc_ld$NB) {
        stop("Error: Length of the weight-at-length vector must equal the number of length bins.")
      }
    }

    if (any(is.na(ma_L))) {
      if (! (is.na(Ls) | is.na(L50)))
        ma_L <- with(blicc_ld, wt_L / (1 + exp(-Ls*(LMP - L50))))    #Mature biomass
      else
        warning("L50 or L95 not specified: mature biomass -at-length not changed.")
    } else {
      if (length(ma_L) != blicc_ld$NB) {
        stop("Error: Length of the mature biomass -at-length vector must equal the number of length bins.")
      }
    }

    if (!is.na(model_name))
      blicc_ld$model_name <- model_name
    if (!is.na(a))
      blicc_ld$a <- a
    if (!is.na(b))
      blicc_ld$b <- b
    if (!is.na(L50))
      blicc_ld$L50 <- L50
    if (!is.na(Ls))
      blicc_ld$Ls <- Ls
    if (!any(is.na(wt_L)))
      blicc_ld$wt_L <- wt_L
    if (!any(is.na(ma_L)))
      blicc_ld$ma_L <- ma_L
    return(blicc_ld)
  }

#' Sets data object from [blicc_dat()] to have a new selectivity priors
#'
#' The priors for each function are defined loosely based on the available
#' data ("empirical Bayes"). The priors are weakly informative, so they set
#' primarily to aid fitting and discourage values outside a reasonable range.
#' For example, the location parameter is generally set with a log-normal
#' mean at the length frequency mode. Estimates far from this point would
#' probably have to be rejected during review after fitting anyway.
#' Prior parameters are selectivity function specific.
#'
#' @inheritParams blicc_mpd
#' @return The data object blicc_ld but with the new selectivity priors
#' @noRd
#'
blip_selectivity <- function(blicc_ld) {
  npar <- Rsel_functions()$npar[blicc_ld$fSel]
  sel_par <- double(sum(npar)) # lognormal mu parameters
  sel_pars <- sel_par          # lognormal sigma parameters

  for (i in 1:blicc_ld$NG) {
    if (blicc_ld$fSel[i]==1) { # logistic
      sm <- with(blicc_ld, c(log((LMP[1]+LMP[which.max(fq[[i]])])*0.75), -1))
    } else if (blicc_ld$fSel[i]==4) {
      sm <- with(blicc_ld, c(log(LMP[which.max(fq[[i]])]), -4, -4))
    } else if (blicc_ld$fSel[i]==5) {
      # This will need reworking...
      sm1 <- with(blicc_ld, LMP[which.max(fq[[i]])])
      if (sm1 > 0.5*(blicc_ld$LMP[1]+blicc_ld$LMP[blicc_ld$NB]))
        sm2 <- with(blicc_ld, 0.5*(LMP[1]+sm1))
      else
        sm2 <- with(blicc_ld, 0.5*(LMP[NB]+sm1))
      sm3 <- 0.2 #with(blicc_ld, 0.02 + fq[[i]][as.integer(sm2)] / max(fq[[i]]))
      sm <- c(log(c(sm1, sm2, sm3)), -5, -5, -5, -5)
      sms <- c(1, 1, 1, 1, 1, 1, 1)
    } else {
      sm <- with(blicc_ld, c(log(LMP[which.max(fq[[i]])]), -4))
    }
    indx <- with(blicc_ld, sp_i[i]:(sp_i[i]+length(sm)-1L))
    sel_par[indx] <- sm
    sel_pars[indx] <- 1.5
  }
  blicc_ld$polSm <- sel_par
  blicc_ld$polSs <- sel_pars
  return(blicc_ld)
}

#' Sets data object from [blicc_dat()] to have a new NB_phi prior
#'
#' The NB_phi lognormal prior is updated with new values mu and sigma.
#' NB_phi controls the negative binomial overdispersion compared to the
#' Poisson distribution, where the variance is mu + mu^2 / NB_phi.
#'
#' @export
#' @inheritParams blip_Linf
#' @param lNBphi  A numeric vector of double containing the mean and sd for the
#'   prior log-normal.
#' @return The data object blicc_ld but with the prior for Galpha changed.
#'
blip_NBphi <- function(blicc_ld,
                       lNBphi) {
  if (!is.numeric(lNBphi) & length(lNBphi)==2)
    stop("Error: lNBphi must be numeric vector of the mu,sd for the NB_phi lognormal prior.")
  blicc_ld$polNB_phim <- lNBphi[1]
  blicc_ld$polNB_phis <- lNBphi[2]
  return(blicc_ld)
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
    nFk    = rep(0.0, blicc_ld$NF),
    nSm    = rep(0, blicc_ld$NP),
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


#' Estimate a minimum number of Gauss-Laguerre quadrature nodes for a defined tolerance
#'
#' The procedure compares population estimates at length to the safe number of nodes
#' to identify the smallest number of nodes that produce estimates that are still
#' within the specified tolerance.
#'
#' @inheritParams blicc_mpd
#' @param draws A dataframe of parameter draws to test
#' @param toler The minimum acceptable tolerance for the integral
#' @return A list of Gauss-Laguerre nodes and weights meeting the minimum tolerance
#' @noRd
#'
LG_Nodes <- function(blicc_ld, draws, toler=1.0e-06) {
  Min_NK <- 50
  glq <-
    statmod::gauss.quad(110, kind = "laguerre", alpha = 0.0)
  suppressWarnings(
    Fk <- as.matrix(tidyr::unnest_wider(dplyr::select(draws, Fk), col="Fk", names_sep="_"))
  )
  suppressWarnings(
    Sm <- as.matrix(tidyr::unnest_wider(dplyr::select(draws, Sm), col="Sm", names_sep="_"))
  )
  Zki <- list()
  pop <- list()
  for (i in 1:nrow(draws)) {
    # Grid values
    Zki[[i]] <- draws$Mk[i] * blicc_ld$M_L
    FSel <- Rselectivities(Sm[i,], blicc_ld)
    for (gi in 1:blicc_ld$NG) {
      if (blicc_ld$Fkg[gi] > 0) {
        Zki[[i]] <- Zki[[i]] + FSel[[gi]] * Fk[i, blicc_ld$Fkg[gi]] # Total mortality
      }
    }
    pop[[i]] <- Cpop_len(glq$nodes, glq$weights, blicc_ld$LLB, Zki[[i]], draws$Galpha[i], draws$Gbeta[i])
  }
  Opt_NK <- Min_NK
  err <- double(1)
  repeat {
    glq <-
      statmod::gauss.quad(Opt_NK, kind = "laguerre", alpha = 0.0)
    for (i in 1:nrow(draws)) {
      err <- max(abs(pop[[i]]-Cpop_len(glq$nodes, glq$weights, blicc_ld$LLB,
                                          Zki[[i]], draws$Galpha[i], draws$Gbeta[i])),
                 na.rm=T)
      if (is.na(err)) stop("Error: Numerical error in LG integration: check the data.")
      else if (err > toler) break
    }
    if (err < toler | Opt_NK >= 110) break
    Opt_NK <- Opt_NK+10
  }
  return(Opt_NK)
}



