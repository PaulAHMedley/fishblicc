# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# OBSERVED-EXPECTED FUNCTIONS
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


# Expected Values ---------------------------------------------------------


#' Generate posterior draws and standard reference points
#'
#' The functions returns a list consisting of the search direction
#' vector, a tibble of draws from the posterior MCMC stanfit (or the
#' mpd estimate) including SPR 20%, 30% and 40% reference points,
#' expected number in each length bin and the data object used to
#' generate these.  The tibble of draws is suitable for use in
#' the `posterior` and `bayesplot` packages and can be used in various
#' functions to show results. Note that the function can take a little time
#' to run and may be moderately large, depending on the number of
#' draws. Some reference points may not exist (indicated by `NA`).
#'
#' @details
#' Reference points are found where possible for the parameter draw from the MCMC.
#' The reference points are found along line representing the relative impact
#' across gears. The default is to apply the same adjustments proportionately
#' to every gear that has a fishing mortality, but alternative scenarios can be
#' defined. For example, if there are two gears, and only the first will be
#' adjusted, then the direction would be `vdir=c(1, 0)`.
#' Where reference points do not exist an `NA` is returned. It is quite
#' possible reference points do not exist in many cases due to the effects of
#' selectivity. Strictly speaking, reference points may not exist within
#' specific constraints on F and selectivity. Within more general bounds on F
#' and selectivity, all reference "points" exist as lines or points. In these
#' functions, reference points in relation to a control (Fk or Sm) are
#' calculated based on the other being fixed at the current estimated value.
#' To see how the fishery might respond when they are adjusted together,
#' the yield or SPR surface must be plotted.
#'
#' @export
#' @inheritParams blicc_mpd
#' @inheritParams blicc_dat
#' @param slimf      An object from the `blicc_fit` or `blicc_mpd` function.
#' @param vdir       A search direction vector with maximum value 1 and
#' minimum 0 applied to changes across gears.
#' @return A list of a tibble containing posterior draws or mpd point estimates
#' of important parameters and per-recruit reference points, the search
#' direction vector and the blicc_ld data object used.
#' @examples
#' \dontrun{
#' res_rp <- blicc_ref_pts(eg_slim, eg_ld)
#' summary(res_rp$rp_df)
#' }
#'
blicc_ref_pts <-
  function(slimf,
           blicc_ld,
           vdir = NA,
           a = NA,
           b = NA,
           L50 = NA,
           L95 = NA) {
    Linf = Galpha = Mk = Fk = Sm = mpd = par = .draw = name2 = se = NULL
    Gbeta = NULL

    vdir <- as.vector(vdir)
    if (is.na(vdir)) {
      vdir <- double(blicc_ld$NG)
      vdir[blicc_ld$Fkg] <- 1
    } else if (max(vdir) != 1.0 | min(vdir) < 0)
      stop(
        "Error: vdir - the direction vector must have at least
            one value equal to 1.0, and no values outside the 0-1 range."
      )

    names(vdir) <- blicc_ld$gname

    if (class(slimf)[1] == "stanfit") {
      # Convert stanfit object to posterior draws with list columns for F and
      # Sm parameters
      rp_df <- posterior::as_draws_df(slimf) |>
        dplyr::select(Linf:`.draw`)

      par_c <- paste0("Sm[", as.character(1:blicc_ld$NP), "]")
      suppressWarnings(
        sel_tmp <- rp_df |>
          dplyr::select(.draw, tidyselect::all_of(par_c)) |>
          tidyr::pivot_longer(
            cols = tidyselect::all_of(par_c),
            names_to = "name2",
            values_to = "Sm"
          ) |>
          dplyr::arrange(.draw, name2) |>
          dplyr::group_by(.draw) |>
          dplyr::summarise(Sm = list(Sm)) |>
          dplyr::ungroup()
      )
      rp_df <- rp_df |>
        dplyr::select(-tidyselect::all_of(par_c))

      par_c <- paste0("Fk[", as.character(1:blicc_ld$NF), "]")
      suppressWarnings(
        F_tmp <- rp_df |>
          dplyr::select(.draw, tidyselect::all_of(par_c)) |>
          tidyr::pivot_longer(
            cols = tidyselect::all_of(par_c),
            names_to = "name2",
            values_to = "Fk"
          ) |>
          dplyr::arrange(.draw, name2) |>
          dplyr::group_by(.draw) |>
          dplyr::summarise(Fk = list(Fk)) |>
          dplyr::ungroup()
      )
      rp_df <- rp_df |>
        dplyr::select(-tidyselect::all_of(par_c))

      rp_df <- rp_df |>
        dplyr::left_join(F_tmp, by = ".draw") |>
        dplyr::left_join(sel_tmp, by = ".draw") |>
        dplyr::select(Linf:Mk, Fk, Sm, tidyselect::everything())

    } else if (class(slimf)[1] == "tbl_df" &&
               all(names(slimf) == c("par", "mpd", "se"))) {
      par_c <- paste0("Sm[", as.character(1:blicc_ld$NP), "]")
      S_v <- dplyr::pull(dplyr::filter(slimf, par %in% par_c), mpd)
      rp_df <- dplyr::filter(slimf,!(par %in% par_c))
      par_c <- paste0("Fk[", as.character(1:blicc_ld$NF), "]")
      F_v <- dplyr::pull(dplyr::filter(rp_df, par %in% par_c), mpd)
      rp_df <- dplyr::filter(rp_df,!(par %in% par_c)) |>
        dplyr::select(-se) |>
        tidyr::pivot_wider(names_from = par, values_from = mpd) |>
        dplyr::mutate(Fk = list(F_v),
                      Sm = list(S_v),
                      `.draw` = 0) |>
        dplyr::select(Linf:Mk, Fk, Sm, tidyselect::everything())

    } else {
      stop("Error: parameter slimf must be a stanfit object or mpd data frame")
    }

    # Update data object
    New_NK <- LG_Nodes(blicc_ld, rp_df)   # Recalculates the LG knots
    if (blicc_ld$NK != New_NK) {
      glq <- statmod::gauss.quad(blicc_ld$NK,
                                 kind = "laguerre", alpha = 0.0)
      blicc_ld$NK <- New_NK
      blicc_ld$gl_nodes <- glq$nodes
      blicc_ld$gl_weights <- glq$weights
    }

    Recalc_SPR <- !all(is.na(c(a, b, L50, L95)))
    if (Recalc_SPR) {
      if (!is.na(a))
        blicc_ld$a <- a
      if (!is.na(b))
        blicc_ld$b <- b
      if (!is.na(L50))
        blicc_ld$L50 <- L50
      if (!is.na(L95))
        blicc_ld$Ls <- -log(1 / 0.95 - 1) / (L95 - blicc_ld$L50)
      blicc_ld$ma_L <- (exp(blicc_ld$b * log(blicc_ld$LMP)) /
                          (1 + exp(
                            -blicc_ld$Ls * (blicc_ld$LMP - blicc_ld$L50)
                          )))
      blicc_ld$wt_L <- blicc_ld$a * exp(blicc_ld$b * log(blicc_ld$LMP))

      rp_df <- rp_df |>
        dplyr::mutate(SPR = purrr::pmap_dbl(
          list(Galpha, Gbeta, Mk, Fk, Sm),
          Calc_SPR,
          blicc_ld = blicc_ld,
          .progress = "SPR"
        ))
    }

    rp_df <- rp_df |>
      dplyr::mutate(
        F20 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          FSPR_solve,
          tarSPR = 0.2,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "F20"
        ),
        F30 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          FSPR_solve,
          tarSPR = 0.3,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "F30"
        ),
        F40 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          FSPR_solve,
          tarSPR = 0.4,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "F40"
        ),
        F01 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          F01_solve,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "F01"
        ),
        S20 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          SSPR_solve,
          tarSPR = 0.2,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "S20"
        ),
        S40 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          SSPR_solve,
          tarSPR = 0.4,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "S40"
        ),
        SMY = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          SMY_solve,
          vdir = vdir,
          blicc_ld = blicc_ld,
          .progress = "SMY"
        )
      )
    lx_df <- blicc_expect_len(rp_df, blicc_ld)
    return(list(
      vdir = vdir,
      rp_df = rp_df,
      lx_df = lx_df,
      ld = blicc_ld
    ))
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
#' @inheritParams blicc_mpd
#' @param rp_df   Posterior draws and reference points tibble from
#' `blicc_ref_pts` function.
#' @return A tibble containing fitted values with respect to length
#' @noRd
#'
blicc_expect_len <- function(rp_df, blicc_ld) {
  Linf = Galpha = Mk = Fk = Sm = .draw = expect = NULL
  suppressWarnings({
    lx_df <- rp_df |>
      dplyr::mutate(expect = purrr::pmap(
        list(Linf, Galpha, Mk, Fk, Sm),
        blicc_get_expected,
        blicc_ld = blicc_ld
      )) |>
      dplyr::select(`.draw`, expect) |>
      tidyr::unnest(expect)
  })
  return(lx_df)
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
#' @inheritParams blicc_mpd
#' @param  Linf     Maximum mean length fish in the population can grow to
#' @param  Galpha   Alpha parameter for the Gamma probability density function
#' that governs growth variability.
#' @param  Mk       Natural mortality (time in units of the growth rate K)
#' @param  Fk       Vector of fishing mortality (time in units of the growth
#' rate K) to be multiplied by selectivity to get fishing mortality in each
#' length bin
#' @param  Sm       Vector of parameters for all the selectivity functions
#' flat-topped selectivity
#' @return A tibble of expected values for each length bin for each gear
#' @examples
#' blicc_get_expected(Linf = 32, Galpha = 100, Mk = 1.5, Fk = 1.5,
#'                    Sm = c(24, 0.1, 0.001), blicc_ld = eg_ld)
#'
blicc_get_expected <-
  function(Linf, Galpha, Mk, Fk, Sm, blicc_ld) {
    # Returns expected values based on the model
    Rsel <- Rselectivities(Sm, blicc_ld)
    pop <- Rpop_F(Galpha, Galpha / Linf, Mk, Fk, Rsel, blicc_ld)

    ex_df <- tibble::tibble()
    for (gi in 1:blicc_ld$NG) {
      efq <- pop$N_L * pop$Fki[[gi]] # Catch
      efq <- sum(blicc_ld$fq[[gi]]) * efq / sum(efq) # Normalise
      ex_df <- rbind(
        ex_df,
        tibble::tibble(
          Sgroup = blicc_ld$gname[gi],
          Lgroup = factor(blicc_ld$LLB),
          sel = Rsel[[gi]],
          N_L = pop$N_L,
          efq = efq
        )
      )
    }
    return(ex_df)
  }


#' A posterior predicted length frequency data set
#'
#' For a set of MCMC parameter draws, returns a simulated catch in each
#' length bin based on the BLICC model. This matrix of simulated data that
#' might be expected if the model is correct, can be used in the posterior
#' and bayesplot packages for evaluation.
#'
#' @export
#' @param blicc_rp List of fishblicc result tables from `blicc_ref_pts()`
#' @param Gear Identifies the single gear (selectivity) providing predictions
#' @param draws  The number of random draws up to the number of
#' draws from the MCMC (the default).
#' @return A matrix with rows equal to gears*draws and columns to length
#' @examples
#' yrep <- posterior_predict(blicc_rp = eg_rp, Gear=1, draws=100)
#'
posterior_predict <- function(blicc_rp, Gear=NA, draws = 0) {
  .draw = Lgroup = Sgroup = efq = NULL
  Gear <- parse_gear(Gear, blicc_rp$ld)
  if (length(Gear) > 1)
    stop("Error: a single gear must be specified.")

  if ((draws <= 0) | (draws >= nrow(blicc_rp$rp_df))) {
    df <- blicc_rp$lx_df
  } else {
    df <- tibble::tibble(`.draw`=sample.int(n = nrow(blicc_rp$rp_df), size = draws)) |>
      dplyr::left_join(blicc_rp$lx_df, by=".draw")
  }

  if (blicc_rp$ld$NG > 1)
    df <- df |>
      dplyr::filter(Sgroup==blicc_rp$ld$gname[Gear])

  # Check order is correct, then extract expected lengths as vector
  ex <- df |>
    dplyr::arrange(.draw, Lgroup) |>
    dplyr::pull(efq)

  phi <- rep(blicc_rp$rp_df$NB_phi, each = blicc_rp$ld$NB)

  yrep <- matrix(
    stats::rnbinom(n = length(ex), mu = ex, size = phi),
    ncol = blicc_rp$ld$NB,
    byrow = TRUE
  )
  # matrix rows will be gears * draws
  return(yrep)
}
