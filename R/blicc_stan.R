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
#' @details
#' Unlike `blicc_fit()`, the model finds a point estimate. This is much faster
#' than the MCMC, but does not provide full information on uncertainty.
#' The model estimates mortality and survival through sequential length
#' intervals. This is then used to derive abundance and catch in numbers within
#' each length bin. The mortality is required to remain constant within each
#' bin, but otherwise can vary arbitrarily. However, in this implementation
#' mortality is constrained by a parametric model with constant natural
#' mortality by length and a double-sided normal selectivity function for
#' fishing mortality. See `blicc_dat()` and `blicc_fit()` documentation for
#' a description of the model and data requirements.
#'
#' The standard errors are estimated from 1000 random draws from a multivariate
#' normal with mean at the MPD mode and covariance estimated from the inverted
#' Hessian matrix. This is an approximation and SE estimates may differ from
#' the full MCMC sampling.
#'
#' In addtion to the parameters, the fit reports the spawning potential ratio
#' (SPR) and the log probability at the mode (lp__).
#'
#' @export
#' @param  blicc_ld    A standard data list suitable for the model
#' (see function `blicc_dat()`)
#' @return A tibble of parameter estimates (mpd) with standard errors.
#' @examples
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0),
#'                 Linf=c(35, 3))
#' mpd_fit <- blicc_mpd(ld)
#'
blicc_mpd <- function(blicc_ld) {
  # Find the posterior mode
  fit <-
    # rstan::optimizing(
    #   stanmodels$BLICC,
    # Test version
    # res <-
    optimizing(
      stmod,
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
  par_names <- c("Linf", "Galpha", "Mk",
                 paste0("Fk", as.character(1:blicc_ld$NF)),
                 paste0("Sm", as.character(1:blicc_ld$NP)),
                 "NB_phi", "Gbeta", "SPR")

  res <- tibble::tibble(
    par = c(full_par_names, "lp__"),
    mpd = unname(c(unlist(fit$par)[par_names],
                   fit$value)),
    se = unname(c(apply(rd, MARGIN=2, FUN="sd"), NA))
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
#' @inheritParams blicc_mpd
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
  # stan_fn <- here("R", "BLICC_mg.stan")
  # stmod <- stan_model(stan_fn, model_name = "BLICC")
  #

  # Find the posterior mode to start
  # # Distributed version
  res <-
    # rstan::optimizing(
    #   stanmodels$BLICC,
    # Test version
    res <-
    optimizing(
      stmod,
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
  # stf <- rstan::sampling(
  #   stanmodels$BLICC,

  # Test version
  stf <- stan(
    fit = stmod,
    file = stan_fn,
    model_name = "BLICC",
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
#' @param LLB  A vector of lower length boundaries for each length frequency
#' bin. Required.
#' @param fq   A list of vectors, each the same length as LLB, containing the
#' frequency data. Zeroes must be included. If only one gear, can be a vector.
#' Required.
#' @param Linf A vector of two values: mean and sd of the prior maximum mean
#' length for the stock. Required.
#' @param sel_fun A vector of selectivity function names (character) or index
#' (1-4), with 1="logistic", 2="normal", 3="ss_normal", 4="ds_normal" for
#' each vector in fq. Required.
#' @param Catch A vector of relative catches. Required, if the number of gears
#' is more than one.
#' @param gear_names A vector of names for the gears for reference. Optional.
#' @param Mk   Natural mortality divided by the growth rate K
#' (usually around 1.5). Optional.
#' @param M_L  Length specific natural mortality vector the same length as
#' Len. Optional.
#' @param a  The length-weight parameter: a*L^b. Not used in model fitting,
#' but for plots etc. Optional.
#' @param b  The length-weight exponent (a*L^b, usually close
#' to 3.0) Optional.
#' @param L50 The length at 50% maturity, often referred to as "first"
#' maturity, and primarily applies to females. Must be less than Linf, usually
#' around 0.66 Linf, which is assumed if it is not provided. Optional.
#' @param L95  The length at 95% maturity. Must be greater than L50 and less
#' than Linf. A small increment is added L50 if it is not provided. Optional.
#' @param sel_par Prior hyper-parameters for the selectivity functions.
#' These are values for the mean of the lognormal for each parameter.
#' Optional.
#' @param spar Matrix of integer indices linking the selectivity function (row)
#' for each parameter in the function (column) to the parameter location.
#' Optional. Must be provided if spar is provided.
#' in the parameter array (e.g. in sel_par). Zeroes mark unused matrix cells.
#' @param NK   Number of nodes for the Gauss Laguerre quadrature rule. 110 is a
#' safe value, but extends the run time of all calculations.  Optional.
#' @return     A list of vectors and numbers structured suitable for use in the
#' Stan model BLICC.stan.
#' @examples
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0), Linf=c(35, 2), sel_fun=4)
#'
blicc_dat <-
  function(LLB,
           fq,
           Linf,
           sel_fun,
           Catch = 1,
           gear_names = NA,
           Mk = NA,
           M_L = 1,
           ma_L = NA,
           wt_L = NA,
           a = 1,
           b = 3,
           L50 = NA,
           L95 = NA,
           sel_par = NA,
           spar = NA,
           NK = NA) {
    # Data vectors
    if (is.vector(fq) & !is.list(fq)) fq <- list(fq)  #convert to list
    if (!is.list(fq)) {
      stop(
        "Error: frequency data must be supplied as a list of vectors (or can be given as a vector if only one.)"
      )
    }
    if (!is.vector(LLB)) {
      stop(
        "Error: please supply the lower bound of bins as a vector of unique vakues and in ascending order."
      )
    }
    NB <- length(LLB)
    LMP <- c((LLB[-NB] + LLB[-1]) * 0.5, LLB[NB] + 0.5*(LLB[NB]-LLB[NB-1])) # Length bin mid points (LMP) used for plotting etc. Not used in Stan model.

    Ngear <- length(fq)
    for (i in 1:Ngear)
      if (!is.vector(fq[[i]]) | (NB != length(fq[[i]]))) {
        stop(
          "Error: please provide the lower bound of bins, and each frequency in each bin, including zeroes as vectors of the same length."
        )
      }

    if (any(diff(LLB) <= 0)) {
      stop("Error: lower bound of bins must be unique and in ascending order.")
    }

    # Selectivity functions
    if (is.character(sel_fun)) {
      sel_fun <- match(sel_fun, c("logistic", "normal", "ss_normal", "ds_normal"))
    }
    if (is.numeric(sel_fun)) {
      if (!(all(sel_fun %in% 1:4))) {
        stop("Error: selectivity functions must be defined as integers or characters: \n 1=logistic, 2=normal, 3=ss_normal, 4=ds_normal.")
      }
    } else {
      stop("Error: selectivity functions must be defined as integers or characters: \n 1=logistic, 2=normal, 3=ss_normal, 4=ds_normal.")
    }
    if (length(sel_fun) != Ngear) {
      stop("Error: the selectivity functions vector (sel_fun) must be the same length as the number of gears.")
    }

    if (is.na(sel_par)) {
      sel_par <- double(Ngear)
      spar <- matrix(0, nrow = Ngear, ncol = 3)
      np <- 1
      for (i in 1:Ngear) {
        if (sel_fun[i] == 4) {
          spar[i,] <- np:(np+2)
          np <- np + 3
        } else {
          spar[i,] <- c(np, np+1, 0)
          np <- np + 2
        }
        if (sel_fun[i]==1) { # logistic
          sm <- c(log((LMP[1]+LMP[which.max(fq[[i]])])*0.75), -1)
        } else if (sel_fun[i]==4) {
          sm <- c(log(LMP[which.max(fq[[i]])]), -4, -4)
        } else {
          sm <- c(log(LMP[which.max(fq[[i]])]), -4)
        }
        indx <- spar[i, spar[i,]>0]
        sel_par[indx] <- sm
      }
      #      sel_par <- sel_prior(LMP, fq, sel_fun, spar)
    } else {
      if (!is.matrix(spar) | any(dim(spar) != c(Ngear, 3)) | length(unique(spar[spar>0])) != max(spar)) {
        stop("Error: log selectivity parameter priors must be indexed in matrix spar. See function details for help.")
      }
      if (length(sel_par) != max(spar) | !is.numeric(sel_par)) {
        stop("Error: log selectivity parameter priors must be defined in vector sel_par and indexed in matrix spar. See function details for help.")
      }
    }
    # Catches
    if (!is.vector(Catch) | length(Catch) != Ngear | any(Catch<0) | !any(Catch>0)) {
      stop("Error: Relative catches must be provided for each length frequency.")
    }
    Fgear <- integer(Ngear)
    gi <- 1
    for (i in 1:Ngear) {
      if (Catch[i] > 0) {
        Fgear[i] <- gi
        gi <- gi + 1
      }
    }

    cprop <- Catch[Catch > 0]/sum(Catch)

    # Growth
    if (length(Linf) < 2) {
      stop("Error: the Linf vector must contain a mean and sd for the prior.")
    }
    if (Linf[1] <= min(LLB)) {
      stop(paste(
        "Error: Linf must be greater than the lowest bin value:",
        as.character(min(LLB))
      ))
    }
    if (Linf[2] <= 0) {
      stop("Error: Linf prior sd must be greater than zero.")
    }
    # Gear names
    if (is.na(gear_names[1])) {
      gear_names <- paste0("Gear_", as.character(1:Ngear))
    }
    if ((length(gear_names) != Ngear) | !is.character(gear_names)) {
      stop("Error: gear_names must be a character vector with the same length as the number of gears.")
    }
    gear_names <- gear_names

    # Weight and maturity
    if (b <= 2 | b > 4) {
      stop("Error: Length-weight exponent (b) must be greater than 2 and less than 4.")
    }
    if (is.na(L50)) {
      L50 <- 0.66 * Linf[1]
    } else {
      if (L50 <= 0.2 * Linf[1] | L50 >= Linf[1]) {
        stop("Error: Length at 50% maturity must be greater than 0.2*Linf and less than Linf.")
      }
    }
    if (is.na(L95)) {
      # Ls = -log(1/0.95 - 1)/(L95-L50)
      Ls <- -log(1 / 0.95 - 1) / (0.05 * (Linf[1] - L50))
    } else {
      if (L95 <= L50 | L95 >= Linf[1]) {
        stop("Error: Length at 95% maturity must be greater than L50 and less than Linf.")
      }
      Ls <- -log(1 / 0.95 - 1) / (L95 - L50)
    }

    if (is.na(ma_L)) {
      ma_L <- (exp(b*log(LMP)) / (1 + exp(-Ls*(LMP - L50))))    #Mature biomass
    } else {
      if (length(ma_L) != NB) {
        stop("Error: Length of the mature biomass -at-length vector must equal the number of length bins.")
      }
    }

    if (is.na(wt_L)) {
      if (a == 1) {
        print("Warning: No weight-at-length information provided - the weight units will be incorrect.")
      }
      wt_L <- a * exp(b*log(LMP))    # Estimated biomass per recruit
    } else {
      if (length(wt_L) != NB) {
        stop("Error: Length of the weight-at-length vector must equal the number of length bins.")
      }
    }

    # Other checks
    if (!is.na(NK))
      if (NK < 50) {
        print(
          "Warning: Having fewer than 50 nodes for the Gauss Laguerre quadarture rule is not advised."
        )
      }

    # Natural mortality
    if (length(M_L) == 1) M_L <- rep(1, NB)
    if (!is.numeric(M_L) | (length(M_L) != NB)) {
      print("Error: The length specific natural mortality vector must be numeric and the same length as length frequencies.")
      return()
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

    dl <- list(
      # Number of knots used in the Gauss-Laguerre quadrature rule
      NK = 110L,
      # Number of gears: separate length frequencies
      NG = Ngear,
      # Number of F's: gears associated with non-zero catches
      NF = length(cprop),
      # Number of length bins
      NB = NB,
      # Number of selectivity parameters
      NP = sum(spar>0),
      # Maximum number of selectivity parameters
      Pmx = 3,
      # Gear names
      gname = as.array(gear_names),
      # Selectivity function indices
      fSel = as.array(sel_fun),
      # Selectivity function parameter index for each gear
      Spar = spar,
      # Estimated total relative catch in numbers of fish, excluding zeroes for surveys etc.
      Catch = as.array(cprop),
      # Index of the Fk associated with each gear: 0 implies catch negligible
      Fkg = as.array(Fgear),
      # Lower length boundaries for each bin
      Len = LLB,
      # Length bin mid points (for plotting etc.)
      LMP = LMP,
      # Frequency data
      fq = fq,
      # Natural mortality variation between lengths vector (Lorenzen M)
      M_L = M_L,
      # Mature biomass at length
      ma_L = ma_L,
      # Biomass at length
      wt_L = wt_L,
      # Length-weight scale parameter
      a = a,
      # Length-weight exponential
      b = b,
      # Length at 50% maturity
      Lm = L50,
      # Slope parameter for the logistic maturity ogive
      Ls = Ls,
      # Expected Linf
      poLinfm = Linf[1],
      # sd for the normal Linf, see above
      poLinfs = Linf[2],
      # Growth mean CV
      polGam = log(1 / 0.1 ^ 2),
      # 10% CV: could be 5% (log(1/0.05^2)) to 30% (log(1/0.3^2)).
      # 30% CV makes length very uninformative however. Default CV=10%
      # growth log-normal sd hyper-parameter
      polGas = 0.25,
      # Natural mortality lognormal mean hyperperamater
      polMkm = lMk[1],
      # Natural mortality lognormal sd hyperperamater
      polMks = lMk[2],
      # Fishing mortality lognormal mean hyperperamater
      # default "fully exploited"
      polFkm = as.array(log(Mk*cprop)),
      # Very weak prior with hyper-parameter lognormal sd
      polFks = 2.0,
      # lognormal mode / 50% selectivity and
      # 1-2 slope parameters for each function
      polSm = sel_par,
      # selectivity lognormal sd
      polSs = rep(1.0, length(sel_par)),
      # Negative binomial lognormal mean
      # Medium level of overdispersion
      polNB_phim = log(100),
      polNB_phis = 0.2,
      #  Catch: sigma for lognormal
      polCs = 0.01
    )
    if (is.na(NK) | NK < 20L) {
      df <- with(dl,
                 expand_grid(
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
      dl$NK <- LG_Nodes(df, dl)
    } else {
      dl$NK <- NK
    }
    glq <-
      statmod::gauss.quad(dl$NK, kind = "laguerre", alpha = 0.0)
    dl$gl_nodes <- glq$nodes
    dl$gl_weights <- glq$weights

    return(dl)
  }


#' Generates selectivity functions priors in a single vector
#'
#' Provides a vector of the prior parameters for selectivity functions based
#' based on the data ("empirical Bayes"). The aim is to provide a loose prior
#' to encourage the MCMC to provide expected values close to the data. In
#' practice the prior should have little influence on the result but support
#' model fitting convergence.
#'
#' @param sel_fun A numeric vector of selectivity functions 1-4
#' @param spar An integer matrix giving the position of each parameter in the return vector
#' @return A vector of the log-parameters for the relevant selectivity
#' functions
#' @noRd
#'
add_sel_prior <- function(blicc_ld) {
  # function for a loose fit of selectivity functions to frequency data
  ls <- function(par, freq, N) {
    ex <- fun(exp(par), blicc_ld$LMP)*RPop
    ex <- ex*N/sum(ex)
    return(sum(((ex-freq)^2)/ex))
  }

  Sv <- RSurvival(gl_nodes, gl_weights, Len, Zki, Galpha, Gbeta)


  res_par <- double(length(fq))
  sm <- double(3)
  for (gi in 1:length(fq)) {
    fun <- switch(sel_fun[gi],
                  Rsel_logistic,
                  Rsel_normal,
                  Rsel_ssnormal,
                  Rsel_dsnormal)

    if (sel_fun[gi]==1) { # logistic
      sm <- c(log((LMP[1]+LMP[which.max(fq[[gi]])])*0.75), 0)
    } else if (sel_fun[gi]==4) {
      sm <- c(log(LMP[which.max(fq[[gi]])]), -4, -4)
    } else {
      sm <- c(log(LMP[which.max(fq[[gi]])]), -4)
    }
    fit <- stats::optim(sm, ls, control = list(maxit=100, reltol = 1e-4), freq=fq[[gi]], N=sum(fq[[gi]]))
    indx <- spar[gi, spar[gi,]>0]
    res_par[indx] <- fit$par
  }
  return(res_par)
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
#' `blicc_ini` if no parameters are provided.
#'
#' @inheritParams blicc_mpd
#' @param  pchain A number indicating the number of chains. Default is 4
#' @param  par    A list of start parameter values (optional)
#' @return A list of start parameters equal to the number of chains for
#' the BLICC Stan model
#' @noRd
#'
blicc_mcmc_ini <- function(blicc_ld, pchain = 4, par = NULL) {
  if (is.null(par)) {
    return(rep(list(blicc_ini(blicc_ld)), pchain))
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
    rstan::extract(stf1,
                   pars = c("lp__"))
  p0 <-
    rstan::extract(stf0,
                   pars = c("lp__"))

  meanp1 <-
    max(p1$lp__) + log(sum(exp(p1$lp__ - max(p1$lp__))) / length(p1$lp__))
  meanp0 <-
    max(p0$lp__) + log(sum(exp(p0$lp__ - max(p0$lp__))) / length(p0$lp__))

  #print("Bayes factor estimate:")
  return(exp(meanp1 - meanp0))
}


#' Estimate a minimum number of Gauss-Laguerre quadrature nodes for a defined tolerance
#'
#' The procedure compares population estimates at length to the safe number of nodes
#' to identify the smallest number of nodes that produce estimates that are still
#' within the specified tolerance.
#'
#' @param draws
#' @param blicc_ld
#' @param toler
#' @return A list of nodes and weights at the minimum level.
#' @NoRd
#'
LG_Nodes <- function(draws, blicc_ld, toler=1.0e-06) {
  Min_NK <- 30
  glq <-
    statmod::gauss.quad(110, kind = "laguerre", alpha = 0.0)
  Fk <- as.matrix(tidyr::unnest_wider(dplyr::select(draws, Fk), col="Fk", names_sep="_"))
  Sm <- as.matrix(tidyr::unnest_wider(dplyr::select(draws, Sm), col="Sm", names_sep="_"))
  Zki <- list()
  pop <- list()
  for (i in 1:nrow(draws)) {
    # Grid values
    Zki[[i]] <- draws$Mk[i] * blicc_ld$M_L
    FSel <- RSelectivities(Sm[i,], blicc_ld)
    for (gi in 1:blicc_ld$NG) {
      if (blicc_ld$Fkg[gi] > 0) {
        Zki[[i]] <- Zki[[i]] + FSel[[gi]] * Fk[i, blicc_ld$Fkg[gi]] # Total mortality
      }
    }
    pop[[i]] <- CPop_Len(glq$nodes, glq$weights, blicc_ld$Len, Zki[[i]], draws$Galpha[i], draws$Gbeta[i])
  }
  Opt_NK <- Min_NK
  err <- double(nrow(draws))
  repeat {
    glq <-
      statmod::gauss.quad(Opt_NK, kind = "laguerre", alpha = 0.0)
    for (i in 1:nrow(draws)) {
      err[i] <- max(abs(pop[[i]]-CPop_Len(glq$nodes, glq$weights, blicc_ld$Len,
                                          Zki[[i]], draws$Galpha[i], draws$Gbeta[i])))
    }
    if (max(err) < toler | Opt_NK >= 110) break;
    Opt_NK <- Opt_NK+10
  }
  return(Opt_NK)
}


#' Returns a tibble containing a summary of the fishblicc priors being applied
#'
#' This is a convenience function that summarizes the priors contained in a
#' fishblicc data object into a tibble. This can be used, for example, to
#' create a table for publication.
#'
#' @export
#' @inheritParams blicc_mpd
#' @return A tibble describing the priors being applied parameter estimates
#' @examples
#' ld <- blicc_dat(LLB = 25:35, fq=c(0,1,2,26,72,66,36,24,12,4,0),
#'                 Linf=c(35, 3))
#' blicc_priors(ld)
#'
blicc_priors <- function(blicc_ld) {

  Selectivity_Functions <- c("Logistic", "Normal", "Single-sided Normal",
                             "Double-sided Normal")
  Sel_par_names <- list(c("Sel 50%", "Sel steepness"), c("Mode", "SD"),
                        c("Mode", "Left SD"), c("Mode", "Left SD", "Right SD"))

  Sel_Func <- Selectivity_Functions[blicc_ld$fSel]

  gear_names <- c(NA, NA, NA, blicc_ld$gname[blicc_ld$Fkg])
  function_type <- c("Normal", rep("Lognormal", 2+blicc_ld$NF))
  sel_par <- vector()
  for (i in 1:blicc_ld$NG) {
    sel_par <- c(sel_par, Sel_par_names[[blicc_ld$fSel[i]]])
    if (dl$fSel[i]==4) {
      gear_names <- c(gear_names, blicc_ld$gname[i], NA, NA)
      function_type <- c(function_type, "Lognormal", NA, NA)
    } else {
      gear_names <- c(gear_names, blicc_ld$gname[i], NA)
      function_type <- c(function_type, "Lognormal", NA)
    }
  }
  gear_names <- c(gear_names, rep(NA, 4))
  function_type <- c(function_type, "Lognormal", rep(NA, 3))

  par_names <- c("Linf", "Galpha", "Mk", rep("Fk", blicc_ld$NF),
                 sel_par, "NB_phi", "b", "Lm", "Ls")

  return(tibble(
    Gear = gear_names,
    Parameter = par_names,
    `Function Type` = function_type,
    Mean = with(blicc_ld, c(poLinfm, exp(c(polGam, polMkm, polFkm, polSm, polNB_phim)),
                            b, Lm, Ls)),
    Mu = with(blicc_ld, c(poLinfm, polGam, polMkm, polFkm, polSm, polNB_phim,
                          b, Lm, Ls)),
    SD = with(blicc_ld, c(poLinfs, polGas, polMks, rep(polFks, NF), polSs, polNB_phis,
                          NA, NA, NA)))
  )

}


