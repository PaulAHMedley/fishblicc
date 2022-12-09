# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# OBSERVED-EXPECTED FUNCTIONS
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


# Expected Values ---------------------------------------------------------


#' Posterior draws including standard reference points for F and Smx
#'
#' The functions returns a tibble of draws from the posterior MCMC stanfit
#' object together with a nested tibble for each draw containing the results
#' (expected values) in terms of length. The tibble is suitable for use in
#' the `posterior` and `bayesplot` packages and can be used in various
#' functions to show results. Note that the function can take a little time
#' to run and may be moderately large, depending on the number of
#' draws, and reference points may not exist (indicated by `NA`).
#'
#' @details
#' Reference points are found where possible for parameter draw from the MCMC.
#' Where reference points do not exist an `NA` is returned. It is quite
#' possible reference points do not exist in many cases due to the effects of
#' selectivity. Strictly speaking, reference points may not exist within
#' specific constraints on F and selectivity. Within more general bounds on F
#' and selectivity, all reference "points" exist as lines or points. In these
#' functions, reference points in relation to a control (Fk or Smx) are
#' calculated based on the other being fixed at the current estimated value.
#' To see how the fishery might respond when they are adjusted together,
#' the yield or SPR surface must be plotted.
#'
#' @export
#' @param slimf      An object from the `blicc_fit` or `blicc_mpd` function.
#' @param blicc_ld   A data list from the `blicc_dat` function.
#' @return A tibble containing posterior draws or mpd point estimates  of
#' important parameters and per-recruit reference points
#' @examples
#' \dontrun{
#' rp_df <- blicc_ref_pts(eg_slim, eg_ld)
#' summary(rp_df)
#' }
#'
blicc_ref_pts <- function(slimf, blicc_ld) {
  Linf=Galpha=Mk=Fk=Smx=Ss1=Ss2=.draw=NULL

  if (class(slimf)=="stanfit") {
    df <- posterior::as_draws_df(slimf)
  } else if (all(class(slimf)==c("tbl_df", "tbl", "data.frame")) &&
             all(mpd$par==c("Linf", "Galpha", "Mk", "Fk", "Smx", "Ss1", "Ss2",
                            "NB_phi", "SPR", "lp__"))) {
    df <- slimf |>
      select(par, mpd) |>
      pivot_wider(names_from=par, values_from=mpd)
  } else {
    return("Error: parameter slimf must be a stanfit object or mpd data frame")
  }

  glq <-
    statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)

  df <- df |>
    dplyr::select(Linf:`.draw`) |>
    dplyr::mutate(
      F20 = purrr::pmap_dbl(
        list(Smx, Linf, Galpha, Mk, Ss1, Ss2),
        FSPR_solve,
        tarSPR = 0.2,
        blicc_ld = blicc_ld,
        glq = glq
      ),
      F30 = purrr::pmap_dbl(
        list(Smx, Linf, Galpha, Mk, Ss1, Ss2),
        FSPR_solve,
        tarSPR = 0.3,
        blicc_ld = blicc_ld,
        glq = glq
      ),
      F40 = purrr::pmap_dbl(
        list(Smx, Linf, Galpha, Mk, Ss1, Ss2),
        FSPR_solve,
        tarSPR = 0.4,
        blicc_ld = blicc_ld,
        glq = glq
      ),
      F01 = purrr::pmap_dbl(
        list(Smx, Linf, Galpha, Mk, Ss1, Ss2),
        F01_solve,
        blicc_ld = blicc_ld,
        glq = glq
      ),
      S20 = purrr::pmap_dbl(
        list(Fk, Linf, Galpha, Mk, Ss1, Ss2),
        SSPR_solve,
        tarSPR = 0.2,
        blicc_ld = blicc_ld,
        glq = glq
      ),
      S40 = purrr::pmap_dbl(
        list(Fk, Linf, Galpha, Mk, Ss1, Ss2),
        SSPR_solve,
        tarSPR = 0.4,
        blicc_ld = blicc_ld,
        glq = glq
      ),
      SMY = purrr::pmap_dbl(
        list(Fk, Linf, Galpha, Mk, Ss1, Ss2),
        SMY_solve,
        blicc_ld = blicc_ld,
        glq = glq
      )
    )
  return(df)
}


#' Generate a data frame of length-based expected values
#'
#' The functions returns a tibble of draws from the posterior MCMC stanfit
#' object together with a nested tibble for each draw containing the results
#' (expected values) in terms of length. The BLICC model and parameter values
#' drawn from an MCMC are used to calculate the expected values for
#' selectivity, survival, relative population numbers, total mortality and
#' expected length frequency for each length bin. The resulting table can be
#' used in various functions to show results. Note that the resulting table
#' of draws may be large depending on the number of draws.
#'
#' @export
#' @param blicc_rp   Posterior draws and reference points tibble from
#' `blicc_ref_pts` function.
#' @param blicc_ld   A data list from the `blicc_dat` function.
#' @return A tibble containing fitted values with respect to length
#' @examples
#' lx_df <- blicc_expect_len(eg_rp, eg_ld)
#' summary(lx_df)
#'
blicc_expect_len <- function(blicc_rp, blicc_ld) {
  Linf=Galpha=Mk=Fk=Smx=Ss1=Ss2=.draw=expect=NULL
  suppressWarnings({
    df <- blicc_rp |>
      dplyr::mutate(expect = purrr::pmap(
        list(Linf, Galpha, Mk, Fk, Smx, Ss1, Ss2),
        blicc_get_expected,
        blicc_ld = blicc_ld
      )) |>
      dplyr::select(`.draw`, expect) |>
      tidyr::unnest(expect)
  })
  return(df)
}


#' Expected length interval selectivity, mortality, survival and frequency
#'
#' For a set of parameter values, returns the expected selectivity, mortality,
#' survival, and catch in each length bin based on the BLICC model as a tibble.
#'
#' @details
#' The model and provided parameter values are used to calculate the expected
#' values for selectivity, survival, relative population numbers,
#' total mortality and expected length frequency for each length bin.
#' Input parameters are point estimates (reals). The standard data list
#' provides the model dimensions and is produced by the function
#' blicc_dat().
#'
#' @export
#' @param  Linf     Maximum mean length fish in the population can grow to
#' @param  Galpha   Alpha parameter for the Gamma probability density function
#' that governs growth variability.
#' @param  Mk       Natural mortality (time in units of the growth rate K)
#' @param  Fk       Fishing mortality (time in units of the growth rate K) to
#' be multiplied by selectivity to get fishing mortality in each length bin
#' @param  Smx      Mode of the normal selectivity function (full selectivity)
#' @param  Ss1      Left side slope, parameterized as 1/sigma^2
#' @param  Ss2      Right side slope, parameterized as 1/sigma^2. Zero implies
#' flat-topped selectivity
#' @param  blicc_ld A standard data list created by `blicc_dat`
#' @return A tibble of expected values for each length bin
#' @examples
#' blicc_get_expected(Linf = 32, Galpha = 100, Mk = 1.5, Fk = 1.5, Smx = 24,
#'                    Ss1 = 0.1, Ss2 = 0.001, blicc_ld = eg_ld)
#'
blicc_get_expected <-
  function(Linf, Galpha, Mk, Fk, Smx, Ss1, Ss2, blicc_ld) {
    # Returns expected values based on the model
    Len <- blicc_ld$Len
    LN <- length(Len)
    glq <-
      statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
    Rsel <- Csel_dsnormal(blicc_ld$LMP, Smx, Ss1, Ss2)
    Fki <- Rsel * Fk
    Zki <- Fki + Mk
    Rsurv <-
      CSurvival_Est(glq$nodes, glq$weights, Len, Zki, Galpha, Galpha / Linf)
    N_L <- CNinInterval(Rsurv, Zki)
    efq <- N_L * Fki # Catch
    efq <- sum(blicc_ld$fq) * efq / sum(efq) # Normalise
    return(tibble::tibble(
      Lgroup = factor(Len),
      sel = Rsel,
      surv = Rsurv,
      N_L = N_L,
      Zk = Zki,
      efq = efq
    ))
  }


#' The expected length frequency
#'
#' For a set of parameter values, returns the expected catch in each length bin
#' based on the BLICC model. This is the same as the blicc_get_expected()
#' function, but only returns the expected length frequency. Used for
#' predictive posterior.
#'
#' @inheritParams blicc_get_expected
#' @param glq  The Gauss-Laguerre nodes and weights: glq <-
#'  statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
#' @return A vector of the expected length bin frequency
#' @examples
#' \dontrun{
#' lfq <- blicc_get_efq(Linf = 32, Galpha = 100, Mk = 1.5, Fk = 1.5, Smx = 24,
#'                      Ss1 = 0.1, Ss2 = 0.001, blicc_ld = eg_ld)
#' plot(y = lfq, x = eg_ld$Len)
#' }
blicc_get_efq <-
  function(Linf, Galpha, Mk, Fk, Smx, Ss1, Ss2, blicc_ld, glq) {
    if (is.na(Fk) | is.na(Smx)) {
      return(NA)
    }
    Len <- blicc_ld$Len
    LN <- length(Len)
    Rsel <- Csel_dsnormal(blicc_ld$LMP, Smx, Ss1, Ss2)
    Fki <- Rsel * Fk
    Zki <- Fki + Mk
    Rsurv <-
      CSurvival_Est(glq$nodes, glq$weights, Len, Zki, Galpha, Galpha / Linf)
    efq <- Fki * CNinInterval(Rsurv, Zki)
    efq <- sum(blicc_ld$fq) * efq / sum(efq) # Normalise
    return(efq)
  }


#' A posterior predicted length frequency data set
#'
#' For a set of MCMC parameter draws, returns a simulated catch in each
#' length bin based on the BLICC model. This matrix of simulated data that
#' might be expected if the model is correct, can be used in the posterior
#' and bayesplot packages for evaluation.
#'
#' @export
#' @inheritParams blicc_expect_len
#' @param draws  The number of random draws up to the number of
#' draws from the MCMC (the default).
#' @return A matrix with rows equal to draws and columns to length
#' @examples
#' yrep <- posterior_predict(blicc_rp = eg_rp, blicc_ld = eg_ld)
#'
posterior_predict <- function(blicc_rp, blicc_ld, draws=0) {
  Linf=Galpha=Mk=Fk=Smx=Ss1=Ss2=NULL
  if ((draws <= 0) | (draws >= nrow(blicc_rp))) {
    df <- blicc_rp
  } else {
    df <- blicc_rp[sample.int(n=nrow(blicc_rp), size=draws),] }

  glq <- statmod::gauss.quad(blicc_ld$NK, kind = "laguerre", alpha = 0.0)
  # Vector of expected frequency at length
  ex <- unlist(purrr::pmap(dplyr::select(df, Linf, Galpha, Mk,
                                         Fk, Smx, Ss1, Ss2),
                           blicc_get_efq, blicc_ld = blicc_ld, glq = glq))
  phi <- rep(df$NB_phi, each = blicc_ld$oBN)

  yrep <- matrix(stats::rnbinom(n=length(ex), mu=ex, size=phi),
                 ncol=blicc_ld$oBN, byrow=TRUE)

  return(yrep)
}
