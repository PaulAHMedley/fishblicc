# BLICC Per Recruit Reference point functions for internal use ---------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

#' Calculate the SPR for the fishing mortalities and selectivities
#'
#' The function calculates the spawning potential ratio (SPR) based on the
#' provided parameter set. The SPR is calculated as a ratio between the spawning
#' biomass per recruit for a particular fishing mortality divided by the
#' spawning biomass per recruit with no fishing. Works with multiple gears and
#' time periods.
#'
#' @inheritParams blicc_get_expected
#' @param Gbeta The Gamma distribution parameter for the growth model
#'   (`Galpha/Linf`)
#' @return The spawning potential ratio
#' @noRd
#' 
Calc_SPR <-
  function(Galpha,
           Gbeta,
           Mk,
           Fk,
           Sm,
           blicc_ld) {
    Sel <- Rselectivities(Sm, blicc_ld)
    SPR <- double(blicc_ld$NT)
    for (ti in seq(blicc_ld$NT)) {
      qi <- blicc_ld$Ti==ti & blicc_ld$Fkq > 0
      selt <- Sel[blicc_ld$Gi[qi]]
      Fkt <- Fk[blicc_ld$Fkq[qi]]
      SPR[ti] <- fSPR(Galpha, Gbeta, Mk, Fkt, selt, blicc_ld)
      }
    SPR0 <- RSPR_0(Galpha, Gbeta, Mk, blicc_ld) # Unexploited SPR
    return(SPR/SPR0)
  }


#' Calculate the YPR for the fishing mortalities and selectivities
#'
#' The function calculates the yield per recruit (YPR) based on the provided
#' parameter set. The YPR is calculated as a sum of catch-at-length multiplied
#' by the weight at length. Works with multiple gears and time periods.
#'
#' @inheritParams Calc_SPR
#' @return The yield per recruit for each gear/time period
#' @noRd
#' 
Calc_YPR <-
  function(Galpha,
           Gbeta,
           Mk,
           Fk,
           Sm,
           blicc_ld) {
    YPR <- double(blicc_ld$NQ)
    Sel <- Rselectivities(Sm, blicc_ld)
    Pop <- Rpop_F(Galpha, Gbeta, Mk, Fk, Sel, blicc_ld)
    for (qi in seq(blicc_ld$NQ)) {
      if (blicc_ld$Fkq[qi] > 0) {
        YPR[qi] <- with(blicc_ld, sum(Pop$N_L[[Ti[qi]]] * Pop$Fki[[Gi[qi]]] * wt_L)) # Catch weight
      }  
    }
    return(YPR)
  }


#' Calculate the relative biomass (depletion)
#'
#' The model calculates the relative biomass compared to the unexploited state
#' by dividing the exploitable biomass in each length bin with
#' fishing by the exploitable biomass with no fishing. The exploitable biomass is
#' the biomass at length weighted by the fishing mortality at each length
#' (i.e. overall selectivity). Works with multiple time 
#' periods.
#'
#' @inheritParams Calc_SPR
#' @return The biomass as a proportion of the unexploited biomass.
#' @noRd
#' 
Calc_BB0 <- function(Galpha, Gbeta, Mk, Fk, Sm, blicc_ld) {
  Bt0 <- double(blicc_ld$NT)
  Sel <- Rselectivities(Sm, blicc_ld)
  popZ <- Rpop_F(Galpha, Gbeta, Mk, Fk,
                 FSel=Sel, blicc_ld)
  Zki <- Mk * blicc_ld$M_L
  popM <- with(blicc_ld,
               Rpop_len(gl_nodes, gl_weights,
                        LLB, Zki, Galpha, Gbeta))
  
  for (ti in seq(blicc_ld$NT)) {
    Fl <- double(blicc_ld$NB)
    for (gi in blicc_ld$Gi[blicc_ld$Ti==ti])
      Fl <- Fl + popZ$Fki[[gi]]
    Wt <- Fl * blicc_ld$wt_L  # Exploitable biomass fish weight * Total F for each length
        
    Bt0[ti] <- sum(popZ$N_L[[ti]] * Wt) / sum(popM * Wt)
  }

  return(Bt0)
}

#' Solve for a fishing mortality which produces the target SPR
#'
#' The function finds the fishing mortality that gives the `tarSPR` spawning
#' potential ratio (SPR). The SPR is calculated as a ratio between the spawning
#' biomass per recruit for a particular fishing mortality divided by the
#' spawning biomass per recruit with no fishing. The method used should
#' generally work for sensible reference point target.
#'
#' However, note that fishing mortality estimate may not be finite, dependent on
#' the selectivity.
#'
#' @details `vdir` should be the same length as gear. It is a weight variable
#'   used to select gears based on values greater than zero. For weights greater
#'   than zero, it implies proportional adjustment to F for each gear relative
#'   to the largest change (1). This is a simple linear change for estimating
#'   reference points. More complex scenarios will need full simulation.
#'
#' @inheritParams blicc_get_expected
#' @param  tarSPR target SPR reference point: usually 0.2, 0.3 or 0.4
#' @param  vdir  A search direction vector with maximum value 1 and minimum 0
#'   applied to changes across gears. Must be the same length as the number of
#'   gear. See details.
#' @return Spawning potential ratio per-recruit fishing mortalities consistent
#'   with parameters and reference point
#' @noRd
#' 
FSPR_solve <-
  function(Linf,
           Galpha,
           Mk,
           Fk,
           Sm,
           tarSPR,
           vdir,
           blicc_ld) {

    SRP_eval <- function(dF) {
      vFk <- (1.0 + vdir*dF)*Fk
      SPR <- fSPR(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
      return(SPR / SPR0 - tarSPR)
    }

    maxval <- (30/max(Fk[which.max(vdir)]) - 1)
    Gbeta <- Galpha / Linf
    Rsel <- Rselectivities(Sm, blicc_ld)
    SPR0 <- RSPR_0(Galpha, Gbeta, Mk, blicc_ld) # Unexploited SPR
    min_SRP <- SRP_eval(-1)
    if (is.na(min_SRP)) return(NA)
    if (min_SRP < 0) return(NA)  # Needs to > tarSPR
    max_SRP <- SRP_eval(maxval)

    if (is.na(max_SRP) | max_SRP > 0) {
      return(NA)
    }
    dF <- stats::uniroot(f = SRP_eval,
                         interval = c(-1, maxval),
                         tol=1e-5,
                         maxiter=500)$root
    return((1.0 + vdir*dF)*Fk)
  }


#' Solve for a selectivity mode which produces the target SPR
#'
#' The function finds the selectivity location parameters that gives the
#' `tarSPR` spawning potential ratio (SPR).
#'
#' Calculates the adjusted location selectivity parameters along a direction
#' vector to achieve a target spawning potential ratio. This should often work
#' for sensible reference point target. However, note that selectivity may not
#' achieve any particular reference point if the fishing mortality is too low.
#' In these cases, `NA` is returned. Only works for a single time period with
#' non-zero F's, so the inputs will need to be filtered accordingly.
#'
#' @details `vdir` should be the same length as gear. It is a dummy variable
#'   used to select gears based on values greater than zero. Unlike for the
#'   fishing mortality, `vdir` is used to select gears where `vdir` > 0. The
#'   value of vdir does not matter. This allows for selectivity mixtures where
#'   changes to selectivity are complex and therefore would need proper
#'   simulation rather than simple linear adjustment available from this
#'   implementation.
#'
#' @inheritParams FSPR_solve
#' @return The selectivity parameter vector with modes (full selectivity)
#'   adjusted to achieve the target SPR. `NA` indicates this target cannot be
#'   achieved.
#' @noRd
#' 
SSPR_solve <-
  function(Linf,
           Galpha,
           Mk,
           Fk,
           Sm,
           tarSPR,
           vdir,
           blicc_ld) {

    SRP_eval <- function(dL) {
      vSm[indx] <- (1 + dL)*Sm[indx]
      Rsel <- Rselectivities(vSm, blicc_ld)
      SPR <- fSPR(Galpha, Gbeta, Mk, Fk, Rsel, blicc_ld)
      return(SPR / SPR0 - tarSPR)
    }

    # select selectivities
    sindx <- get_selectivities(which(vdir>0), blicc_ld)
    indx <- blicc_ld$sp_i[sindx]

    ref_par <- max(Sm[indx]) # location par
    maxdL <- Linf/ref_par - 1

    Gbeta <- Galpha / Linf
    vSm <- Sm
    SPR0 <- RSPR_0(Galpha, Gbeta, Mk, blicc_ld) # Unexploited SPR

    # First need to bracket tarSPR
    S2 <- maxdL
    V2 <- SRP_eval(S2)
    if (V2 < 0) {
      S1 <- maxdL
      V1 <- V2
      S2 <-
        Linf * (1 + 5 / sqrt(Galpha)) # A length that fish do not grow to
      S2 <- S2/ref_par - 1
      V2 <- SRP_eval(S2)
      if (V2 < 0) {
        return(NA)
      }
    } else {
      # V2 > 0
      mindL <- blicc_ld$LLB[1]/ref_par - 1
      S1 <- blicc_ld$L50/ref_par - 1
      V1 <- SRP_eval(S1)
      while ((V1 > 0) & (S1 >= mindL)) {
        S1 <- S1 - 1.0
        V1 <- SRP_eval(S1)
      }
      if (V1 > 0) {
        return(NA)
      }
    }
    dL <- stats::uniroot(f = SRP_eval,
                         interval = c(S1, S2),
                         tol=1e-5,         # Depends on binwidth precision
                         maxiter=500)$root
    vSm[indx] <- (1 + dL)*Sm[indx]
    return(vSm)
  }


#' Find the standard yield-per-recruit fishing mortality reference point F0.1
#'
#' Estimates an F0.1 reference point consistent with parameters and YPR slope
#' 10% of initial slope. The function solves for F0.1 using a simple linear
#' approximation to slope at the origin and at each point on the curve. Only
#' works for a single time period with non-zero F's, so the inputs will need to
#' be filtered accordingly.
#'
#' @inheritParams FSPR_solve
#' @return F0.1 per-recruit fishing mortality reference point
#' @noRd
#' 
F01_solve <-
  function(Linf,
           Galpha,
           Mk,
           Fk,
           Sm,
           vdir,
           blicc_ld) {
    fYPRF01 <- function(dF) {
      vFk <- (1.0 + vdir*dF)*Fk
      vFkd <- (1.0 + vdir*(dF+0.01))*Fk
      # Crude approx: slope=(y2-y1)/(x2-x1). constant (x2-x1) so ignore
      slope <-
        (
          fYPR(Galpha, Gbeta, Mk, vFkd, Rsel, blicc_ld) -
            fYPR(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
        )
      return(slope - slope01)
    }

    maxval <- (30/max(Fk[which.max(vdir)]) - 1)
    Gbeta <- Galpha / Linf
    Rsel <- Rselectivities(Sm, blicc_ld)

    vFk <- (1.0 - vdir)*Fk
    vFkd <- (1.0 - vdir*0.99)*Fk
    slope01 <-
      0.1 * (
        fYPR(Galpha, Gbeta, Mk, vFkd, Rsel, blicc_ld) -
          fYPR(Galpha, Gbeta, Mk, vFk, Rsel, blicc_ld)
      )
    if (fYPRF01(maxval) > 0) {
      return(NA)
    }
    dF <- stats::uniroot(f = fYPRF01,
                         interval = c(-1, maxval),
                         tol=1e-5,
                         maxiter=500)$root
    return((1.0 + vdir*dF)*Fk)
  }

#' Solve for a selectivity mode which produces maximum yield
#'
#' The function finds the selectivity mode (`Smx`) that maximises the yield per
#' recruit (YPR), keeping all other parameters at their current value.
#'
#' This should generally work because there must be a maximum yield between
#' the extreme lengths (`0` and `Linf`). Only works for a single 
#' time period with non-zero F's, so the inputs will need to be filtered 
#' accordingly. 
#'
#' @inheritParams FSPR_solve
#' @return maximum yield point for the selectivity mode
#' @noRd
#'
SMY_solve <-
  function(Linf, Galpha, Mk, Fk, Sm, vdir, blicc_ld) {

    YPR_S <- function(dL) {
      vSm[indx] <- (1 + dL)*Sm[indx]
      Rsel <- Rselectivities(vSm, blicc_ld)
      Yield <- fYPR(Galpha, Gbeta, Mk, Fk, Rsel, blicc_ld)
      return(Yield)
    }

    Gbeta <- Galpha / Linf
    vSm <- Sm
    # index location parameters for gears contributing to fishing mortality
    # select selectivities
    sindx <- get_selectivities(which(vdir>0), blicc_ld)
    indx <- blicc_ld$sp_i[sindx]

    maxdL <- Linf/max(Sm[indx]) - 1

    res <- stats::optimize(
      YPR_S,
      lower = -1,
      upper = maxdL,
      maximum = TRUE,
      tol = 1.0e-5)

    vSm[indx] <- (1 + res$maximum)*Sm[indx]
    return(vSm)
  }


#' Yield per recruit with variable fishing mortality. 
#'
#' Used to evaluate the YPR for different values of Fk. Only works for a single 
#' time period with non-zero F's, so the inputs will need to be filtered 
#' accordingly. 
#'
#' @inheritParams FSPR_solve
#' @param Gbeta  Rate parameter for the Gamma distribution growth variability
#'   (`=Galpha/Linf`)
#' @param RSel   List of selectivities
#' @return yield-per-recruit
#' @noRd
#'
fYPR <- function(Galpha, Gbeta, Mk, Fk, Rsel, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L
  Fki <- list()
  for (gi in seq(blicc_ld$NG)) {
    Fki[[gi]] <- Rsel[[gi]] * Fk[gi]
    Zki <- Zki + Fki[[gi]]      # Total mortality
  }

  N_L <- with(blicc_ld, Cpop_len(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta))
  Yield <- 0
  for (gi in seq(blicc_ld$NG)) 
    Yield <- Yield + sum(N_L * Fki[[gi]] * blicc_ld$wt_L)
  return(Yield)
}


#' Spawning biomass per recruit with variable fishing mortality
#'
#' Used to evaluate the SPR for different values of Fk.
#'
#' @inheritParams fYPR
#' @return spawning biomass per recruit
#' @noRd
#'
fSPR <- function(Galpha, Gbeta, Mk, Fk, Rsel, blicc_ld) {
  # Spawning potential ratio
  Zki <- Mk * blicc_ld$M_L
  for (gi in seq_along(Fk)) {
    Zki <- Zki + Rsel[[gi]] * Fk[gi]      # Total mortality
  }
  N_L <- with(blicc_ld, Cpop_len(gl_nodes, gl_weights, LLB,
                                 Zki, Galpha, Gbeta))
  return(sum(N_L * blicc_ld$ma_L))
}


#' Optimum length with maximum yield
#'
#' Used to find the maximum YPR (at optimum exploitation length) subject to a
#' SPR constraint. NOT USED
#'
#' @inheritParams fYPR
#' @param tarSPR  Target Spawner Potential Ratio
#' @return maximum yield-per-recruit
#' @noRd
#' 
maxYPR <- function(Galpha, Gbeta, Mk, tarSPR, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L

  N_L <- with(blicc_ld, Cpop_len(gl_nodes, gl_weights, LLB, Zki, Galpha, Gbeta))
  Pwt_L <- N_L  * blicc_ld$wt_L
  return(sum(Pwt_L))
}

