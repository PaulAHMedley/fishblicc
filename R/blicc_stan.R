# BLICC Stan fitting functions



blicc_fit <- function(ld,
                      ptargn = 2000,
                      pwup = 1000,
                      pthin = 1,
                      pchain = 4) {
  #' Fits a Bayesian length interval catch curve to a length frequency data sample
  #'
  #' Uses a constant recruitment model with double-sided normal selectivity, a single natural
  #' mortality parameter, and Gamma-distributed variable growth to define expected numbers of fish
  #' in pre-defined length bins. The model assumes a negative binomial likelihood function to fit
  #' to a single length frequency sample. Selectivity, mortality, asymptotic mean length and scale
  #' parameters are fitted. In addition, the spawning potential ratio is generated which ca be used
  #' to indicate the resource status.
  #'
  #' @export
  #' @param  ld A standard data list (see function BLICC_dat)
  #' @param  ptargn target draw for the MCMC
  #' @param  pwup   warm up iterations for the Stan MCMC
  #' @param  pthin  thinning for the MCMC draws - only set above 1 if there is significant autocorrelation
  #' @param  pchain number of chains to run in parallel
  #' @return An object of class `stanfit` returned by `rstan::sampling`

  # @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
  # Find the posterior mode to start
  res <-
    rstan::optimizing(
      stanmodels$BLICC,
      data = ld,
      init = blicc_ini(),
      hessian = TRUE,
      as_vector = FALSE,
      verbose = FALSE,
      iter = 20000,
      refresh = 500,
      tol_obj = 1e-12,
      tol_rel_obj = 1e3,
      tol_grad = 1e-8,
      tol_rel_grad = 1e5,
      tol_param = 1e-8
    )

  niter <- pwup + pthin * ptargn / pchain
  stf <- rstan::sampling(stanmodels$BLICC,
    data = ld,
    chains = pchain,
    control = list(adapt_delta = 0.90, max_treedepth = 12),
    iter = niter,
    warmup = pwup,
    thin = pthin,
    init = blicc_mcmc_ini(pchain, res$par),
    verbose = FALSE
  )
  return(stf)
}  # blicc_fit


blicc_dat <-
  function(LLB,
           fq,
           Linf,
           lMk = NA,
           b = 3,
           L50 = 0.66,
           Ls = 0.05) {
    #' Generates a data list, including first-guess suitable parameter values, for the BLICC model
    #'
    #' Vectors for the lower bound for each length bin and the number of fish in each length bin and combined with
    #' prior parameters into a data list as expected by the BLICC model. Constant bin widths are
    #' not necessary for the model and not assumed. Zero frequencies must be explicitly in the fq vector.
    #'
    #' @export
    #' @param LLB   A vector of the lower length boundaries for the length frequency
    #' @param fq    A vector of the same length as LLB containing the frequency data. Zeroes must be included.
    #' @param Linf  A vector of two values being the mean and sd of the prior mean maximum length for the species
    #' @param lMk   Log natural mortality divided by the growth rate K (usually around 1.5)
    #' @param b     A number with the length-weight exponent (usually close to 3.0)
    #' @param L50   A number for the length at 50% maturity as a proportion of Linf (primarily applies to females). Must be less than 1.0.
    #' @param Ls    A number for the logistic steepness parameter for the maturity ogive
    #' @param bw    A number providing the bin width (usually 1)
    #' @return      A list of vectors and numbers structured suitable for use in the Stan model BLICC.stan
    #'

    # check input vectors, are the same length, and LLB must be in sequence and parameters are in range
    if (!(is.vector(LLB) &
          is.vector(fq)) | (length(LLB) != length(fq))) {
      print(
        "Error: please supply two vectors of same length 1: lower bound of bins 2: frequency in each bin"
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
    if (b <= 2) {
      print("Error: Length-weight exponent (b) must be greater than 2.")
      return()
    }
    if (L50 <= 0.2 | L50 >= 1.0) {
      print("Error: Length at 50% maturity must be greater than 0.2 and less than 1.0.")
      return()
    }
    if (Ls <= 0.0) {
      print("Error: Maturity ogive slope (Ls) must be greater than 0.")
      return()
    }
    if (bw <= 0) {
      print("Error: bin width must be greater than zero.")
      return()
    }

    # First guess for maximum selectivity
    Smx = c((lfd$LLB[which(max(lfd$fq) == lfd$fq)])[1], 5)                     # For double-sided normal
    #  S50 = c(((lfd$LLB[1] + lfd$LLB[which(max(lfd$fq) == lfd$fq)])[1])*0.5, 5)  # For logistic

    if (is.na(lMk)) {
      lMk <- c(log(b * (1 - L50) / L50),   #from Prince et al. 2015
               0.1)
    }

    return(
      list(
        NK = 110,
        oBN = length(lfd$fq),
        Len = lfd$LLB,
        fq = lfd$fq,
        b = b,
        Lm = L50 * Linf[1],
        Ls = Ls * Linf[1],
        # The Linf prior
        poLinfm = Linf[1],
        poLinfs = Linf[2],
        polGam = log(1 / 0.1 ^ 2),
        polGas = 0.25,
        polMkm = lMk[1],
        polMks = lMk[2],
        polFkm = log(1.5),
        polFks = 2.0,
        polSmxm = log(Smx[1]),
        polSmxs = log(Smx[2]),
        polSs1m = -0.1,
        polSs1s = 2.0,
        polSs2m = -5.0,
        polSs2s = 2.0,
        polNB_phim = log(100),
        polNB_phis = 0.5
      )
    )
  } # blicc_dat



blicc_ini <- function() {
  #' Default start parameters centred on priors
  #'
  #' @return A list of initial parameter values for the BLICC Stan model
  #'
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
} # blicc_ini


blicc_mcmc_ini <- function(pchain=4, par = NULL) {
  #' Expands start parameters to the number of chains required for the MCMC run
  #'
  #' @param  pchain A number indicating the number of chains. Default is 4
  #' @param  par    A list of start parameter values (optional)
  #' @return A list of start parameters equal to the number of chains for the BLICC Stan model
  #'
  if (is.null(par)) {
    return(rep(list(BLITCC_Ini()), pchain))
  } else {
    return(rep(list(par), pchain))
  }
} # blicc_mcmc_ini

