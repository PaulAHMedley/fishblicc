# BLICC Model Functions ---------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><

#' Returns a list of supported selectivity functions with parameters
#'
#' The function documents the list of available selectivity
#' functions. This is used in various routines as a
#' constant.
#'
#' @return  A list of function names, short names and parameter list.
#' @noRd
#'
Rsel_functions <- function() {
  sel <- list(long_name = c("Logistic",
                      "Normal",
                      "Single-sided Normal",
                      "Double-sided Normal"),
        short_name = c("logistic",
                       "normal",
                       "ssnormal",
                       "dsnormal"),
        par_names  = list(
                       logistic = c("Sel 50%", "Steepness"),
                       normal = c("Mode", "SD"),
                       ssnormal = c("Mode", "Left SD"),
                       dsnormal = c("Mode", "Left SD", "Right SD"))
        )
   sel$npar <- sapply(sel$par_names, FUN=length)
   return(sel)
}


#' Parse the selectivity function, converting to an integer if necessary
#'
#' The selectivity function parameter is converted to an integer index of the function, or `NA`
#' is returned with an error message.
#'
#' @inheritParams blicc_mpd
#' @param sel_fun  A number or string representing a valid available function
#' @return Integer gear index or NA if no gear can be specified
#' @noRd
#'
parse_selectivity <- function(sel_fun, blicc_ld) {
  func <- Rsel_functions()
  Nfunc <- length(func$short_name)
  if (is.character(sel_fun)){
    if (! all(sel_fun %in% func$short_name))
      stop(paste("Error: specified selectivity function must be the function name (",
                  func$short_name,
                  ") or an integer between 1 and ", as.character(Nfunc)))
    sel_fun <- match(sel_fun, func$short_name)
  } else {
    if (! is.numeric(sel_fun))
      stop(paste0("Error: specified selectivity function must be the function name (",
                  func$short_name,
                  ") or an integer between 1 and", as.character(Nfunc)))
    sel_fun <- as.integer(sel_fun)
    if (min(sel_fun) < 1 | max(sel_fun) > Nfunc)
      stop(paste("Error: specified selectivity function must be the function name (",
                  func$short_name,
                  ") or an integer between 1 and", as.character(Nfunc)))
  }
  return(sel_fun)
}


#' Parse the gear parameter, converting to an integer
#'
#' The gear parameter is converted to an integer index of the gear, or `NA`
#' is returned with an error message.
#'
#' @inheritParams blicc_mpd
#' @param Gear  A number or string representing a valid gear in the model
#' @return Integer gear index or NA if no gear can be specified
#' @noRd
#'
parse_gear <- function(Gear, blicc_ld) {
  if (blicc_ld$NG == 1) {
    return(1)
  } else {
    if (any(is.na(Gear))) {
      stop(
        paste0(
          "Error: Gears must be specified that exactly matches a gear name or be an integer between 1 and ",
          as.character(blicc_ld$NG)
        )
      )
    } else {
      if (is.character(Gear[1])) {
        Gear <- match(Gear, blicc_ld$gname)
        if (any(is.na(Gear))) {
          stop(
            paste0(
              "Error: Specified gears must exactly match a gear name or be an integer between 1 and ",
              as.character(blicc_ld$NG)
            )
          )
        }
      }
      Gear <- as.integer(Gear)
      if (!all(dplyr::between(Gear, 1, blicc_ld$NG))) {
        stop(paste0(
          "Error: Specified gear must be between 1 and ",
          as.character(blicc_ld$NG)
        ))
      }
      return(Gear)
    }
  }
}


# Model 1
#' Calculate a logistic selectivity curve for a length vector
#'
#' A logistic takes a 50% selectivity and a slope parameter. The resulting
#' selectivity varies from 0 to 1.0, with 1.0 being the asymptote. It is
#' calculated for each position defined in the length vector.
#'
#' @export
#' @param  Sp   Vector (length 2) of the selectivity parameters (50%,
#'              and slope) for the logistic selectivity function.
#' @param  LMP  A vector of lengths (usually mid-points for length bins)
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_logistic(Sp = c(30, 0.5), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_logistic <- function(Sp, LMP) {
  SL <- double(length(LMP))
  diff <- LMP - Sp[1L]
  SL <- 1.0/(1.0 + exp(-Sp[2L] * diff))
  return(SL)
}

# Model 2
#' Calculate a normal selectivity curve for a length vector
#'
#' A normal takes one slope (1/sigma^2) parameter for the
#' both sides of the mode parameter (mu). The resulting selectivity varies
#' from 0 to 1.0, with 1.0 being the mode. It is calculated for each position
#' defined in the length vector.
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 2) of the selectivity parameters (mode,
#'              and 1/sigma^2) for the normal selectivity function.
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_normal(Sp = c(35, 0.2), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_normal <- function(Sp, LMP) {
  SL <- double(length(LMP))
  diff <- (LMP - Sp[1L]) ^ 2
  SL <- exp(-Sp[2L] * diff)
  return(SL)
}

# Model 3
#' Calculate a single-sided normal selectivity curve for a length vector
#'
#' A single-sided normal takes one slope (1/sigma^2) parameter for the
#' left side of the mode parameter (mu). For the right side the slope is assumed
#' to be zero. The resulting selectivity varies from 0 to 1.0,
#' with 1.0 being the asymptote. It is calculated for each position defined
#' in the length vector.
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 2) of the selectivity parameters (mode,
#'              and left side 1/sigma^2) for the normal selectivity function.
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_ssnormal(Sp = c(35, 0.2), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_ssnormal <- function(Sp, LMP) {
  SL <- double(length(LMP))
  S1 <- LMP < Sp[1L]
  SL[S1] <- exp(-Sp[2L] * ((LMP[S1] - Sp[1L]) ^ 2))
  SL[!S1] <- 1.0
  return(SL)
}

# Model 4
#' Calculate a double-sided normal selectivity curve for a length vector
#'
#' A double-sided normal takes two separate slope (1/sigma^2) parameters
#' around the mode parameter (mu) that change the slope independently
#' on either side of the mode. The selectivity varies from 0 to 1.0,
#' with 1.0 being at the mode. It is calculated for each position defined
#' in the length vector. Setting the second slope parameter to zero
#' creates a flat-top selectivity (see `Rsel_ssnormal`).
#'
#' @export
#' @inheritParams Rsel_logistic
#' @param  Sp   Vector (length 3) of the selectivity parameters (mode,
#'              left  and right side 1/sigma^2) for the normal
#'              selectivity function
#' @return      A vector of selectivity values varying from 0.0 to 1.0
#' @examples
#' Sel <- Rsel_dsnormal(Sp = c(35, 0.1, 0.2), LMP = seq(15.5, 55.5, by=1.0))
#' plot(y=Sel, x=seq(15.5, 55.5, by=1.0))
#'
Rsel_dsnormal <- function(Sp, LMP) {
  # Double sided normal
  SL <- double(length(LMP))
  S1 <- LMP < Sp[1L]
  SL[S1] <- exp(-Sp[2L] * (LMP[S1] - Sp[1L]) ^ 2)
  SL[!S1] <- exp(-Sp[3L] * (LMP[!S1] - Sp[1L]) ^ 2)
  return(SL)
}


#' Calculate the survival of a fish cohort to sequential length boundaries
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the
#' survival to each length bin boundary taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. Survival can be used to estimate population
#' numbers-at-length and catch length frequency, among other things.
#' The model uses the Gauss-Laguerre quadrature rule for integration,
#' so nodes and weights for this must be provided.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param Galpha Alpha parameter for the Gamma probability density function
#' that governs growth variability.
#' @param Gbeta  Rate parameter for the Gamma distribution growth variability
#' (=Galpha/Linf)
#' @param Mk     Natural mortality divided by the growth rate K
#' @param Fk     Fishing mortality divided by the growth rate K for each gear
#' making a contribution
#' @param FSel    A list of all the selectivities for each length bin
#' @return A list of the population size in each length bin and a list of
#' fishing mortalities at length for each gear.
#' @examples
#' Sel <- Rselectivities(exp(eg_ld$polSm), eg_ld)
#' S <- Rpop_F(100, 100/50, Mk=1.5, Fk=exp(eg_ld$polFkm),
#'             FSel=Sel, blicc_ld=eg_ld)
#' plot(y=S$NL, x=eg_ld$LMP, type="l")
#'
Rpop_F <- function(Galpha, Gbeta, Mk, Fk, FSel, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L
  for (gi in 1:blicc_ld$NG) {
    if (blicc_ld$Fkg[gi] > 0) {
      FSel[[gi]] <- FSel[[gi]] * Fk[blicc_ld$Fkg[gi]]  # Fishing mortality
      Zki <- Zki + FSel[[gi]]                          # Total mortality
    }
  }

  N_L <- with(blicc_ld,
              Cpop_len(gl_nodes, gl_weights,
                       LLB, Zki, Galpha, Gbeta) )
  return(list(N_L=N_L, Fki=FSel))
}


#' Calculate the survival of a fish cohort to sequential length boundaries
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the
#' relative population number in each length bin taking into account
#' variation in growth and varying mortality-at-length. This is
#' analogous to an age based catch curve model, but for length.
#' The model uses the Gauss-Laguerre quadrature rule for integration,
#' so nodes and weights for this must be provided.
#'
#' @export
#' @inheritParams Rpop_F
#' @param  node   Nodes for the Gauss-Laguerre quadrature rule
#' @param  wt     Weights for the Gauss-Laguerre quadrature rule
#' @param  Len    Vector of lower-bounds for length intervals
#' @param  Zki    Vector of total mortality (in time units of growth rate
#' parameter K) to be applied in each length interval.
#' This vector should be the same length as Len.
#' @return A vector of proportion surviving to each length bin lower bound.
#' @examples
#' glq <- statmod::gauss.quad(90, kind = "laguerre", alpha = 0.0)
#' S <- Rpop_len(glq$nodes, glq$weights, Len=15:55,
#'               Zki=c(rep(1.5, 10), rep(3, 31)), 100, 100/50)
#' plot(y=S, x=15:55, type="l")
Rpop_len <- function(node, wt, Len, Zki, Galpha, Gbeta)  {
  zsum <- function(x, Lr, Z) {
    # used to apply sequenced sum of mortality
    return(sum(log(x + Lr) * Z))
  }
  lgamma_Galpha <- lgamma(Galpha)
  nv <- length(node)
  LN <-  length(Len)
  surv <- double(LN)
  x_beta <- node / Gbeta
  log_x_beta <- log(x_beta)

  ss <- log(node + Gbeta * Len[1]) * (Galpha - 1.0) -
    Gbeta * Len[1] - lgamma_Galpha
  surv[1] <- sum(exp(ss) * wt)
  ss <-
    (log(x_beta + Len[2] - Len[1]) * (-Zki[1])) +
    log_x_beta * Zki[1] +
    log(node + Gbeta * Len[2]) * (Galpha - 1.0) -
    Gbeta * Len[2] - lgamma_Galpha
  surv[2] <- sum(exp(ss) * wt)

  for (Li in 3:LN) {
    Ln <- Len[Li]
    Lrange <- Ln - Len[1:(Li - 1)]
    Zii <- c(-Zki[1], Zki[1:(Li - 2)] - Zki[2:(Li - 1)])
    v2 <- log_x_beta * Zki[Li - 1]
    lim <- vapply(
      x_beta,
      zsum,
      FUN.VALUE = 0,
      Lr = Lrange,
      Z = Zii,
      USE.NAMES = FALSE
    )
    ss <- lim + v2 + log(node + Gbeta * Ln) * (Galpha - 1.0) -
      Gbeta * Ln - lgamma_Galpha
    surv[Li] <- sum(exp(ss) * wt)
  }
  pop <- c(surv[-LN] - surv[-1], surv[LN]) / Zki
  return(pop)
}


#' Calculate the survival of a fish cohort to sequential length boundaries
#'
#' The model integrates over growth rate variability and sums mortality
#' piece-wise over intervening length intervals to derive the
#' survival to each length bin boundary taking into account variation in growth
#' and varying mortality-at-length. This is analogous to an age based catch
#' curve model, but for length. Survival can be used to estimate population
#' numbers-at-length and catch length frequency, among other things.
#' The model uses the Gauss-Laguerre quadrature rule for integration,
#' so nodes and weights for this must be provided.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param Sm  Vector of selectivity parameters for all the
#' selectivity functions
#' @return A list of fishing mortalities at length for each gear.
#' @examples
#' S <- Rselectivities(exp(eg_ld$polSm), eg_ld)
#' plot(y=S[[1]], x=eg_ld$LMP, type="l")
Rselectivities <- function(Sm, blicc_ld) {
  Ski <- list()
  for (gi in 1:blicc_ld$NG) {
    Indx <- blicc_ld$spar[gi,]
    Indx <- Indx[Indx>0]
    #
    # Need to be changed to C functions or a C function
    #
    Ski[[gi]] <- switch(blicc_ld$fSel[gi],
                        Rsel_logistic(Sm[Indx], blicc_ld$LMP),
                        Rsel_normal(Sm[Indx], blicc_ld$LMP),
                        Rsel_ssnormal(Sm[Indx], blicc_ld$LMP),
                        Rsel_dsnormal(Sm[Indx], blicc_ld$LMP))
  }
  return(Ski)
}

#' Calculate the spawning potential with no fishing
#'
#' The model calculates the spawning potential by multiplying the proportion of
#' numbers in each length bin with no fishing by the mature biomass per recruit.
#'
#' @export
#' @inheritParams Rpop_F
#' @return The spawning potential (a double) when F=0.
#' @examples
#' RSPR_0(Galpha=100, Gbeta=100/50, Mk=1.5, blicc_ld=eg_ld)
#'
RSPR_0 <- function(Galpha, Gbeta, Mk, blicc_ld) {
  Zki <- Mk * blicc_ld$M_L
  pop <- with(blicc_ld,
              Rpop_len(gl_nodes, gl_weights,
                       LLB, Zki, Galpha, Gbeta))
  return(sum(pop * blicc_ld$ma_L))
}


