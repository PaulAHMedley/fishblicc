# BLICC Observed-Expected Functions -------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

#' Generate posterior draws and standard reference points
#'
#' The function returns a list consisting of the search direction vector, a
#' tibble of draws from the posterior MCMC stanfit (or the mpd estimate)
#' including SPR 20% and 40% reference points, expected number in each length
#' bin and the data object used to generate these. The tibble of draws is
#' suitable for use in the posterior and bayesplot packages and can be used in
#' various functions to show results. Note that the function can take a little
#' time to run and may be moderately large, depending on the number of draws.
#' Some reference points may not exist (indicated by `NA`).
#'
#' @details The reference points calculated are F40 and F20 for 40% and 20% SPR
#'   respectively, and changes in the selectivity location parameter to achieve
#'   either SPR 40% (S40) or maximum yield based on yield-per-recruit (SMY). As
#'   well as reference points, the function recalculates the SPR if the
#'   length-weight or maturity parameters are provided. In addition, the
#'   exploitable biomass as a proportion of the unexploited biomass (B_B0) for
#'   each time period and the YPR for each data frequency are calculated for
#'   each parameter draw.
#'
#'   Reference points are calculated for a particular scenario, represented by a
#'   time period. By default, this time period is the last in the data set,
#'   although the period can be specified. This is necessary because different
#'   time periods can have different gears operating, which affects fishing
#'   mortality at length and therefore the fishing mortality reference points.
#'   This has no effect if there is only one time period.
#'
#'   Reference points are found where possible for the parameter draw from the
#'   MCMC. The reference points are found along a line representing the relative
#'   impact across gears. The default is to apply the same adjustments
#'   proportionately to every gear that has a fishing mortality, but alternative
#'   scenarios can be defined. For example, for two gears the default is
#'   `vdir=c(1, 1)`. If there are two gears, and only the first will be
#'   adjusted, then the direction would be `vdir=c(1, 0)`. It is possible to
#'   provide a vector such as `vdir=c(1, 0.5)` to represent only a 50% change in
#'   the second gear compared to the first, for example. This will only be
#'   applied to the fishing mortality. The selectivity changes will continue to
#'   use the vector as a dummy variable, but the quantity will be ignored and
#'   all `vdir` with non-zero values for gears will be treated as having equal
#'   displacement; so for example a selectivity evaluation for `vdir = c(1.0,
#'   0.5)` is equivalent to `vdir = c(1, 1)`.
#'
#'   It is likely that proper simulation projections would be required for more
#'   complex cases where gears are not being simply excluded from management
#'   action. Unless there is a specific need, the best option is to leave it as
#'   the default `vdir=NULL` so the reference points are indicative and
#'   consistent.
#'
#'   Where reference points do not exist an `NA` is returned. It is quite
#'   possible reference points do not exist in many cases due to the effects of
#'   selectivity. Strictly speaking, reference points may not exist within
#'   specific constraints on F and selectivity. Within more general bounds on F
#'   and selectivity, all reference "points" exist as lines or points. In these
#'   functions, reference points in relation to a control (`Fk` or `Sm`) are
#'   calculated based on the other being fixed at the current estimated value.
#'   To see how the fishery might respond when they are adjusted together, the
#'   yield or SPR surface must be plotted.
#'
#' @export
#' @inheritParams blicc_mpd
#' @inheritParams blicc_dat
#' @param slimf An object from the [blicc_fit] or [blicc_mpd] functions.
#' @param time_period Specifies the time period used for reference point
#'   scenario. If not specified, the default is the final time period.
#' @param vdir  A search direction vector with maximum value 1 and minimum 0
#'   applied to changes across gears. Optional.
#' @return A list of 1) dr_df: posterior draws or mpd point estimates of
#'   important parameters 2) lx_df: observed-expected length frequencies, 3) ld:
#'   `blicc_ld` data list used. 4) scenario used: includes a direction vector
#'   and subsets of the data and parameters sufficient to calculate rp_df 5)
#'   rp_df: per-recruit reference points
#' @examples
#' \dontrun{
#' res_rp <- blicc_ref_pts(trgl_slim, trgl_ld)
#' summary(res_rp$rp_df)
#' }
#' 
blicc_ref_pts <-
  function(slimf,
           blicc_ld,
           time_period = NULL,
           vdir = NULL,
           a = NA,
           b = NA,
           L50 = NA,
           L95 = NA) {
    Linf = Galpha = Mk = Fk = Sm = mpd = par = .draw = name2 = se = NULL
    Gbeta = SPR = lp__ = B_B0 = NULL

    if (is.null(time_period)) {
      time_period <- blicc_ld$NT
    } else {
      time_Period <- parse_period(time_period, blicc_ld)
      if (length(time_period) > 1) 
        stop("Error: Only one time period can be specified to define the reference point scenario.\n")
    }
    
    gears <- with(blicc_ld, Gi[Ti==time_period & Fkq > 0])
    
    if (length(gears)==0) 
      stop("Error: A time period must be specified with at least one active fishing gear.\n")
    
    if (is.null(vdir)) {
      vdir <- rep(1, length(gears))
    } else if (length(vdir) != length(gears) | max(vdir) != 1.0 | min(vdir) < 0)
      stop(
        paste("Error: vdir - the direction vector must have at least
            one value equal to 1.0, and no values outside the 0-1 range and be of length", 
              as.character(length(gears)), "\n")
      )

    names(vdir) <- blicc_ld$gname[gears]

    if (class(slimf)[1] == "stanfit") {
      # Convert stanfit object to posterior draws with list columns for F and
      # Sm parameters
      dr_df <- posterior::as_draws_df(slimf) |>
        dplyr::select(Linf:`.draw`)

      par_c <- paste0("Sm[", as.character(1:(blicc_ld$NP+blicc_ld$NM)), "]")
      suppressWarnings(
        sel_tmp <- dr_df |>
          dplyr::select(.draw, tidyselect::all_of(par_c)) |>
          tidyr::pivot_longer(
            cols = tidyselect::all_of(par_c),
            names_to = "name2",
            values_to = "Sm"
          ) |>
          dplyr::group_by(.draw) |>
          dplyr::summarise(Sm = list(Sm)) |>
          dplyr::ungroup()
      )
      dr_df <- dr_df |>
        dplyr::select(-tidyselect::all_of(par_c))

      par_c <- paste0("Fk[", as.character(1:blicc_ld$NF), "]")
      suppressWarnings(
        F_tmp <- dr_df |>
          dplyr::select(.draw, tidyselect::all_of(par_c)) |>
          tidyr::pivot_longer(
            cols = tidyselect::all_of(par_c),
            names_to = "name2",
            values_to = "Fk"
          ) |>
          dplyr::group_by(.draw) |>
          dplyr::summarise(Fk = list(Fk)) |>
          dplyr::ungroup()
      )
      dr_df <- dr_df |>
        dplyr::select(-tidyselect::all_of(par_c))

      par_c <- paste0("SPR[", as.character(1:blicc_ld$NT), "]")
      suppressWarnings(
        SPR_tmp <- dr_df |>
          dplyr::select(.draw, tidyselect::all_of(par_c)) |>
          tidyr::pivot_longer(
            cols = tidyselect::all_of(par_c),
            names_to = "name2",
            values_to = "SPR"
          ) |>
          dplyr::group_by(.draw) |>
          dplyr::summarise(SPR = list(SPR)) |>
          dplyr::ungroup()
      )
      dr_df <- dr_df |>
        dplyr::select(-tidyselect::all_of(par_c))
      
      dr_df <- dr_df |>
        dplyr::left_join(F_tmp, by = ".draw") |>
        dplyr::left_join(sel_tmp, by = ".draw") |>
        dplyr::left_join(SPR_tmp, by = ".draw") |>
        dplyr::select(Linf:Mk, Fk, Sm, tidyselect::everything())

    } else if (class(slimf)[1] == "tbl_df" &&
               all(names(slimf) == c("par", "mpd", "se"))) {
      par_c <- paste0("Sm[", as.character(1:(blicc_ld$NP+blicc_ld$NM)), "]")
      S_v <- dplyr::pull(dplyr::filter(slimf, par %in% par_c), mpd)
      dr_df <- dplyr::filter(slimf,!(par %in% par_c))
      par_c <- paste0("Fk[", as.character(1:blicc_ld$NF), "]")
      F_v <- dplyr::pull(dplyr::filter(dr_df, par %in% par_c), mpd)
      dr_df <- dplyr::filter(dr_df,!(par %in% par_c))
      par_c <- paste0("SPR[", as.character(1:blicc_ld$NT), "]")
      R_v <- dplyr::pull(dplyr::filter(dr_df, par %in% par_c), mpd)
      
      dr_df <- dplyr::filter(dr_df,!(par %in% par_c)) |>
        dplyr::select(-se) |>
        tidyr::pivot_wider(names_from = par, values_from = mpd) |>
        dplyr::mutate(Fk = list(F_v),
                      Sm = list(S_v),
                      SPR = list(R_v),
                      `.draw` = 0) |>
        dplyr::select(Linf:Mk, Fk, Sm, tidyselect::everything())
    } else {
      stop("Error: parameter slimf must be a stanfit object or mpd data frame")
    }

    # Update data object
    New_NK <- LG_Nodes(blicc_ld, dr_df)   # Recalculates the LG knots
    if (blicc_ld$NK != New_NK) {
      glq <- statmod::gauss.quad(New_NK,
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

      dr_df <- dr_df |>
        dplyr::mutate(SPR = purrr::pmap(
          list(Galpha, Gbeta, Mk, Fk, Sm),
          Calc_SPR,
          blicc_ld = blicc_ld,
          .progress = "SPR"
        ))
    }

    dr_df <- dr_df |>
      dplyr::mutate(
        B_B0 = purrr::pmap(
          list(Galpha, Gbeta, Mk, Fk, Sm),
          Calc_BB0,
          blicc_ld = blicc_ld,
          .progress = "B_B0"
        ),
        YPR = purrr::pmap(
          list(Galpha, Gbeta, Mk, Fk, Sm),
          Calc_YPR,
          blicc_ld = blicc_ld,
          .progress = "YPR"
      )) |>
      dplyr::select(Linf:lp__, SPR, B_B0, YPR, dplyr::everything())
    
    lx_df <- blicc_expect_len(dr_df, blicc_ld)
    
    # For reference points, we only estimate for the reference time_period scenario    
    Findx <- blicc_ld$Ti==time_period & blicc_ld$Fkq > 0
    seli <- get_selectivities(blicc_ld$Gi[Findx], blicc_ld)
    parindx <- integer(0)
    for (si in seli) {
      parindx <- with(blicc_ld, c(parindx, sp_i[si]:sp_e[si]))
    }
    # Add mixtures
    for (gi in blicc_ld$Gi[Findx]) {
      gii <- gi*2L - 1L
      if (blicc_ld$GSmix1[gii] > 0)
        parindx <- with(blicc_ld, c(parindx, NP + GSmix1[gii]:GSmix1[gii+1L]))
    }

    # subset dr_df for the tp_ld data
    rp_df <- dr_df |>
      dplyr::select(-lp__, -SPR, -B_B0, -YPR) |>
      dplyr::mutate(
        Fk = purrr::pmap(
          list(Fk),
          \(x) x[Findx]
        ),
        Sm = purrr::pmap(
          list(Sm),
          \(x) x[parindx]
        )
      )

    tp_ld <- blicc_period_filter(blicc_ld, time_period) # reduces ld to single time period
    
    rp_df <- rp_df |>
      dplyr::mutate(
        F20 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          FSPR_solve,
          tarSPR = 0.2,
          vdir = vdir,
          blicc_ld = tp_ld,
          .progress = "F20"
        ),
        # F30 = purrr::pmap(
        #   list(Linf, Galpha, Mk, Fk, Sm),
        #   FSPR_solve,
        #   tarSPR = 0.3,
        #   vdir = vdir,
        #   blicc_ld = blicc_ld,
        #   .progress = "F30"
        # ),
        F40 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          FSPR_solve,
          tarSPR = 0.4,
          vdir = vdir,
          blicc_ld = tp_ld,
          .progress = "F40"
        ),
        # F01 = purrr::pmap(
        #   list(Linf, Galpha, Mk, Fk, Sm),
        #   F01_solve,
        #   vdir = vdir,
        #   blicc_ld = blicc_ld,
        #   .progress = "F01"
        # ),
        # S20 = purrr::pmap(
        #   list(Linf, Galpha, Mk, Fk, Sm),
        #   SSPR_solve,
        #   tarSPR = 0.2,
        #   vdir = vdir,
        #   blicc_ld = blicc_ld,
        #   .progress = "S20"
        # ),
        S40 = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          SSPR_solve,
          tarSPR = 0.4,
          vdir = vdir,
          blicc_ld = tp_ld,
          .progress = "S40"
        ),
        SMY = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          SMY_solve,
          vdir = vdir,
          blicc_ld = tp_ld,
          .progress = "SMY"
        )
      )
    
    return(list(
      dr_df = dr_df,
      lx_df = lx_df,
      ld = blicc_ld,
      scenario = list(time_period = time_period,
                      gears = gears,
                      vdir = vdir,
                      time_period_ld = tp_ld
                      ),
      rp_df = rp_df
    ))
  }


#' Generate expected proportional catch numbers for each gear.
#'
#' The function returns a tibble of the expected relative catch number for each
#' gear for the provided parameters. The gear values add up to 1, so for single
#' gear fisheries, this is will return 1; therefore this is only useful in
#' multi-gear fisheries.
#'
#' As well as the expected catch, the fit residuals are returned as a
#' diagnostic. For the mpd fit, the function returns a tibble with the number of
#' rows equal to the number of gears. For the MCMC fit, it returns the number of
#' MCMC iterations multiplied by the number of gears.
#'
#' @export
#' @param blicc_rp  A list of posterior draws, reference points with associated
#'   direction, the data object and expected lengths from [blicc_ref_pts]
#'   function.
#' @return A tibble of expected catches and residuals for each gear and MCMC
#'   iteration
#' @examples
#' blicc_expected_catches(trgl_rp) |>
#'   dplyr::group_by(gear) |>
#'   dplyr::summarise(resid=mean(resid))
#' 
blicc_expected_catches <- function(blicc_rp) {
  Linf = Galpha = Mk = Fk = Sm = .draw = NULL
  
  dr_df <- blicc_rp$dr_df
  blicc_ld <- blicc_rp$ld  
  suppressWarnings({
    ca_df <- dr_df |>
      dplyr::mutate(
        expect = purrr::pmap(
          list(Linf, Galpha, Mk, Fk, Sm),
          blicc_get_eca,
          blicc_ld = blicc_ld),
        resid = purrr::pmap(
          list(expect),
          \(x) as.vector(blicc_ld$prop_catch - x))
      ) |>
      dplyr::select(`.draw`, expect, resid) |>
      tidyr::unnest_longer(c(expect, resid)) |> 
      dplyr::mutate(
        gear = rep(blicc_ld$gname, nrow(dr_df)),
        std_resid = resid / blicc_ld$polCs) |>
      dplyr::select(`.draw`, gear, expect, resid, std_resid) 
  })
  return(ca_df)
}


#' Generate a data frame of length-based expected values
#'
#' The functions returns a tibble of draws from the posterior MCMC stanfit
#' object together with a nested tibble for each draw containing the results
#' (expected values) in terms of length. The BLICC model and parameter values
#' drawn from an MCMC are used to calculate the expected values for selectivity,
#' survival, relative population numbers, total mortality and expected length
#' frequency for each length bin. The resulting table can be used in various
#' functions to show results. Note that the resulting data frame may be
#' large depending on the number of draws.
#'
#' @inheritParams blicc_mpd
#' @param dr_df   Posterior draws and reference points tibble from
#'   [blicc_ref_pts] function.
#' @return A tibble containing fitted values with respect to length
#' @noRd
#'
blicc_expect_len <- function(dr_df, blicc_ld) {
  Linf = Galpha = Mk = Fk = Sm = .draw = expect = NULL
  suppressWarnings({
    lx_df <- dr_df |>
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


#' Get expected BLICC selectivity, mortality, survival and length frequency
#'
#' For a set of parameter values, returns the expected selectivity, mortality,
#' survival, and catch in each length bin based on the BLICC model as a tibble.
#'
#' @details The model and provided parameter values are used to calculate the
#'   expected values for selectivity, survival, relative population numbers,
#'   total mortality and expected length frequency for each length bin. Input
#'   parameters are point estimates. A standard data list provides the
#'   model dimensions as produced by the function [blicc_dat].
#'
#' @export
#' @inheritParams blicc_mpd
#' @inheritParams Rselectivities
#' @param  Linf     Maximum mean length fish in the population can grow to
#' @param  Galpha   Alpha parameter for the Gamma probability density function
#'   that governs growth variability.
#' @param  Mk       Natural mortality (time in units of the growth rate K)
#' @param  Fk       Vector of fishing mortality (time in units of the growth
#'   rate K) to be multiplied by selectivity to get fishing mortality in each
#'   length bin
#' @return A tibble of expected values for each length bin for each gear
#' @examples
#' blicc_get_expected(Linf = 62, Galpha = 100, Mk = 1.5, Fk = 1.5,
#'                    Sm = c(24, 0.1), blicc_ld = gillnet_ld)
#'
blicc_get_expected <-
  function(Linf, Galpha, Mk, Fk, Sm, blicc_ld) {
    # Returns expected values based on the model
    Rsel <- Rselectivities(Sm, blicc_ld)
    pop <- Rpop_F(Galpha, Galpha / Linf, Mk, Fk, Rsel, blicc_ld)

    ex_df <- tibble::tibble()
    for (qi in seq(blicc_ld$NQ)) {
      gi <- blicc_ld$Gi[qi]
      ti <- blicc_ld$Ti[qi]
      efq <- pop$N_L[[ti]] * pop$Fki[[gi]] # Catch
      efq <- sum(blicc_ld$fq[[qi]]) * efq / sum(efq) # Normalise
      ex_df <- with(blicc_ld, rbind(
        ex_df,
        tibble::tibble(
          Qgroup = fqname[qi],
          Lgroup = factor(LLB),
          sel = Rsel[[gi]],
          N_L = pop$N_L[[ti]],
          efq = efq
        )
      ))
    }
    return(ex_df)
  }


#' Get the expected length frequency from BLICC model parameters 
#'
#' For a set of parameter values, returns the expected catch in each length bin
#' for each gear based on the BLICC model. This is the same as the
#' [blicc_get_expected] function, but only returns the expected length
#' frequency. Used for predictive posterior in some plots. 
#'
#' @inheritParams blicc_get_expected
#' @param gear_i A integer vector of gears to obtain the expected frequency for
#' @return A list of two vectors: the gear names, and a vector of the expected
#'   length bin frequency
#' @noRd
#'
blicc_get_efq <-
  function(Linf, Galpha, Mk, Fk, Sm, gear_i, blicc_ld) {
    if (any(is.na(Fk)) | any(is.na(Sm))) {
      return(NA)
    }
    Rsel <- Rselectivities(Sm, blicc_ld)
    Pop <- Rpop_F(Galpha, Galpha/Linf, Mk, Fk, Rsel, blicc_ld)

    efq <- double(0)
    for (qi in seq(blicc_ld$NQ)) {
      gi <- blicc_ld$Gi[qi]
      if (gi %in% gear_i) {
        ti <- blicc_ld$Ti[qi]
        ex_fq <- Pop$N_L[[blicc_ld$Ti[qi]]] * Pop$Fki[[gi]] # Catch
        ex_fq <- sum(blicc_ld$fq[[gi]]) * ex_fq / sum(ex_fq) # Normalise
        efq <- c(efq, ex_fq)
      }
    }
    # interleaved gear_names=rep(blicc_ld$gear_names[gear_i], each=blicc_ld$NB)
    return(efq)
  }


#' Get the expected catch proportions from BLICC model parameters
#'
#' For a set of parameter values, returns the expected catch numbers for each
#' gear for each gear based on the BLICC model. This is the same as the
#' [blicc_get_expected] function, but only returns the expected catches. 
#'
#' @inheritParams blicc_get_expected
#' @return A vector of the same length as the number of frequencies containing
#'   the expected catch number proportions.
#' @noRd
#' 
blicc_get_eca <-
  function(Linf, Galpha, Mk, Fk, Sm, blicc_ld) {
    if (blicc_ld$NG == 1) {
      return(rep(1, blicc_ld$NQ))
    }
    Rsel <- Rselectivities(Sm, blicc_ld)
    pop <- Rpop_F(Galpha, Galpha/Linf, Mk, Fk, Rsel, blicc_ld)
    
    eca <- double(blicc_ld$NQ)
    for (qi in seq(blicc_ld$NQ)) {
      eca[qi] <- with(blicc_ld, sum(pop$N_L[[Ti[qi]]] * pop$Fki[[Gi[qi]]])) # Catch numbers
    }
    sums <- with(blicc_ld, tapply(eca, Ti[Fkq>0], sum))
    eca <- with(blicc_ld, eca / sums[Ti[Fkq>0]]) # Normalise
    return(as.vector(eca))
  }


#' Get a posterior predicted length frequency
#'
#' For a set of MCMC parameter draws, returns a simulated catch in each
#' length bin based on the BLICC model. This matrix of simulated data that
#' might be expected if the model is correct, can be used in the posterior
#' and `bayesplot` packages for evaluation.
#'
#' @export
#' @param blicc_rp List of fishblicc result tables from [blicc_ref_pts]
#' @param gear Identifies the single gear (selectivity) providing predictions
#' @param time_period Identifies the single time_period used to simulate a 
#'   length frequency
#' @param draws  The number of random draws up to the number of draws from the
#'   MCMC (the default).
#' @return A matrix with rows equal to gears*draws and columns to length
#' @examples
#' yrep <- posterior_predict(blicc_rp = trgl_rp, gear=1, time_period = 1, draws=100)
#'
posterior_predict <- function(blicc_rp, gear = NULL, time_period = NULL, draws = 0) {
  .draw = Lgroup = Qgroup = efq = NULL
  
  blicc_ld <- blicc_rp$ld
  gear <- parse_gear(gear, blicc_ld)
  if (length(gear) > 1)
    stop("Error: a single gear must be specified. \n")
  time_period <- parse_period(time_period, blicc_ld)
  if (length(time_period) > 1)
    stop("Error: a single time_period must be specified. \n")
  
  fqi <- with(blicc_ld, which(gear==Gi & time_period==Ti))
  
  if ((draws <= 0) | (draws >= nrow(blicc_rp$dr_df))) {
    df <- blicc_rp$lx_df
  } else {
    df <- tibble::tibble(`.draw`= sample.int(n = nrow(blicc_rp$dr_df), 
                                             size = draws)) |>
      dplyr::left_join(blicc_rp$lx_df, by=".draw")
  }

  if (blicc_rp$ld$NQ > 1)
    df <- df |>
      dplyr::filter(Qgroup==blicc_ld$fqname[fqi])

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
