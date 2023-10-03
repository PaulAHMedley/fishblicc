# BLICC Functions for Setting Priors -----------------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><



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
      warning(paste0("Mk outside range 1-5: ",
                     format(exp(lMk[1]), digits=2), " (make sure you provide the log-Mk)"))
    blicc_ld$polMkm <- lMk[1]
  }
  if (!is.na(lMk[2]))
    blicc_ld$polMks <- lMk[2]
  return(blicc_ld)
}


#' Sets data object from [blicc_dat] to have a new Fk prior
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
  if (any(is.na(lFks)))
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


#'Sets data object to have a new life-history vectors
#'
#'Sets the maturity at length and weight-at-length vectors in the data object
#'(from [blicc_dat]) to have new values. The values are either provided as
#'vectors or model parameters are provided to calculate the vectors.
#'Calculations are done for the length bin mid-points. The models for the
#'calculation are the standard length-weight (W=aL^b)  and logistic model for
#'the maturity-at-length (Lm = 1/(1+exp(Ls*(L-L50)))).
#'
#'@export
#'@inheritParams blip_Linf
#'@inheritParams blicc_dat
#'@param set_defaults Logical indicating whether to set defaults or not if
#'  parameters are NA. Leaves then alone if FALSE.
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
        warning("No weight-at-length information provided - the weight units will be incorrect.")
        a <- 1.0
      }
      if (is.na(b) & set_defaults)
        b <- 3.0
      else if (b <= 2 | b > 4) {
        stop("Error: Length-weight exponent (b) must be greater than 2 and less than 4.")
      }
      if (!(is.na(a) | is.na(b)))
        wt_L <-
          with(blicc_ld, a * exp(b * log(LMP)))    # Estimated biomass per recruit
      else
        warning("a or b not specified: weight-at-length not changed.")
    } else {
      if (length(wt_L) != blicc_ld$NB) {
        stop("Error: Length of the weight-at-length vector must equal the number of length bins.")
      }
    }
    
    if (any(is.na(ma_L))) {
      if (!(is.na(Ls) | is.na(L50)))
        ma_L <-
          with(blicc_ld, wt_L / (1 + exp(-Ls * (LMP - L50))))    #Mature biomass
      else
        warning("L50 or L95 not specified: mature biomass -at-length not changed.")
    } else {
      if (length(ma_L) != blicc_ld$NB) {
        stop(
          "Error: Length of the mature biomass -at-length vector must equal the number of length bins."
        )
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

#' Sets data object from [blicc_dat] to have a new selectivity priors
#'
#' The priors for each function are defined loosely based on the available data
#' ("empirical Bayes"). It is recommended that priors are weakly informative, so
#' they are set primarily to aid fitting and discourage values outside a
#' reasonable range.
#'
#' The prior hyper-parameters are "estimated" if there is a single selectivity
#' for a gear. The estimates are based on the 50% (median), 10% and 90%
#' quartiles of the cumulative frequency. This ensures that the selectivity
#' function is centred on the frequency data and provides a robust estimate of
#' reasonable priors. Alternatively, setting parameters manually can be done
#' using the [blip_set_sel] function.
#'
#' If a selectivity is made up of a mixture of functions, hyper-parameters are
#' not estimated and [blip_set_sel] must be used to set the hyperparameters
#' manually. This is necessary because the user must propose the hypothesis for
#' the mixtures and it is too easy for the selectivity model to be overfitted
#' otherwise.
#'
#' @export
#' @inheritParams blicc_selfun
#' @return The data object blicc_ld but with the new selectivity priors
blip_selectivity <- function(blicc_ld,
                             sel_indx = NULL) {
  if (is.null(sel_indx))
    sel_indx = 1L:blicc_ld$NS
  else 
    sel_indx <- parse_sel_indx(sel_indx, blicc_ld)

  # ssd <- 1.281552
  ssd <- 2.0 # sd parameter for selectivity slopes 10%-90% range
  Galpha <- exp(blicc_ld$polGam)
  Linf <- blicc_ld$poLinfm
  Gbeta <- Galpha/Linf
  LLB <- blicc_ld$LLB
  LMP <- blicc_ld$LMP
  Zk <- exp(blicc_ld$polMkm)*blicc_ld$M_L
  gl <- statmod::gauss.quad(110, "laguerre", alpha=0)
  pop <- Rpop_len(gl$nodes, gl$weights, LLB, Zk, Galpha, Gbeta)

  for (si in sel_indx) {
    par_range <- with(blicc_ld, sp_i[si]:sp_e[si])
    #    if (any(is.na(blicc_ld$polSm[par_range]))) {  # Check any previous estimate to preserve
    # find this selectivity function in a gear as a single function
    ggi <- which(blicc_ld$GSbase == si)
    gi <- NA
    if (length(ggi) > 0) {
      ggi <- ggi[blicc_ld$GSmix1[2L * (ggi - 1L) + 1L] == 0]
      gi <- ggi[1]
    }
    
    if (is.na(gi))
      warning(
        paste0(
          "Selectivity function ",
          as.character(si),
          " is in a mixture, so the prior will need to be set directly. ",
          "See function `blip_set_sel`. \n"
        )
      )
    else {
      pfq <- blicc_ld$fq[[gi]] / pop  # adjust data for mortality
      pfq <- pfq / sum(pfq)           # normalise
      cfq <- cumsum(pfq)            # cumulative sum
      i10 <-
        max(1L, which(cfq < 0.1), na.rm = TRUE)    # 10% quartile index
      i25 <-
        min(which(cfq > 0.25), blicc_ld$BN - 1L, na.rm = TRUE)   # 25% quartile index
      i50 <-
        min(which(cfq > 0.5), blicc_ld$BN - 1L, na.rm = TRUE)    # 50% quartile index
      i90 <-
        min(which(cfq > 0.9), blicc_ld$BN - 1L, na.rm = TRUE)    # 90% quartile index
      
      switch(blicc_ld$fSel[si],
             {
               #logistic
               blicc_ld$polSm[par_range] <-
                 log(c(LMP[i25],
                       abs(0.5 * (
                         log(exp(0.1) - 1) / (LMP[i10] - LMP[i25]) +
                           log(exp(0.9) - 1) / (LMP[i90] - LMP[i25])
                       ))))
               # slopes very loosely based on integral of the logistic
             },
             {
               #normal
               blicc_ld$polSm[par_range] <-
                 log(c(LMP[i50], (0.5 * (LMP[i90] - LMP[i10]) / ssd) ^
                         -2))
             },
             {
               #ssnormal
               blicc_ld$polSm[par_range] <-
                 log(c(LMP[i50], ((LMP[i50] - LMP[i10]) / ssd) ^ -2))
             },
             {
               #dsnormal
               blicc_ld$polSm[par_range] <-
                 log(c(LMP[i50], 
                       ((LMP[i50] - LMP[i10]) / ssd) ^ -2, 
                       ((LMP[i90] - LMP[i50]) / ssd) ^ -2))
             })
      blicc_ld$polSs[par_range] <- 1.5  # default
    }
  }
  #  }
  return(blicc_ld)
}


#' Sets a selectivity function's prior hyper-parameters to particular values
#'
#' Each parameter set is for a single selectivity function indexed by
#' `sel_indx`. Although the hyper-parameters can be set directly in the data
#' object, this function does the indexing making it a little easier. The
#' location parameter is set as a length measure. The slope parameters are set
#' as log values.
#'
#' The location parameter (`loc`) is either the 50% selectivity for the logistic
#' or the mode of the particular normal selectivity function. The `lslope`
#' parameter is a single log value of the slope for the relevant function. This
#' is either the log of the logistic slope parameter, or the log of the
#' reciprocal of the variance for the normal functions. For the double-sided
#' normal, two parameters must be provided and for all other functions only one
#' parameter. If `lslope` is not provided, a default slope of -4.5 is applied if
#' the slope parameters have not already set. If values are already present
#' these are conserved.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param sel_indx Single integer indexing a selectivity function
#' @param loc Single location parameter for the selectivity function
#' @param lslope Vector of prior selectivity log slope parameters manually set for
#'   the model. Only one gear at a time can be specified to use this. Parameters
#'   must be in the correct order: in function order, with each function having
#'   a location parameter followed by one or two steepness parameters.
#' @return The data object `blicc_ld` but with the new selectivity priors set
#'
blip_set_sel <- function(blicc_ld,
                         sel_indx,
                         loc,
                         lslope = NULL) {
  
  sel_indx <- parse_sel_indx(sel_indx, blicc_ld, TRUE)

  if (! (is.vector(loc) & is.numeric(loc)) |
      length(loc) != 1 |
      loc <= 0)
    stop("Error: `loc` must be a single positive numeric value for logistic 50% selectivity or normal mode. \n")

  if (blicc_ld$fSel[sel_indx] == 4) np <- 2 else np <- 1

  if (! is.null(lslope)) {
    if (! (is.vector(lslope) & is.numeric(lslope)) |
        length(lslope) != np)
      stop(paste0("Error: `lslope` for this selectivity function must be a numeric vector of length ",
                  as.character(np), " being the log of the slope parameter. \n"))
  }

  blicc_ld$polSm[blicc_ld$sp_i[sel_indx]] <- log(loc)
  if (is.null(lslope))
    blicc_ld$polSm[(blicc_ld$sp_i[sel_indx]+1):blicc_ld$sp_e[sel_indx]] <- -4.5 #default
  else
    blicc_ld$polSm[(blicc_ld$sp_i[sel_indx]+1):blicc_ld$sp_e[sel_indx]] <- lslope
  if (any(is.na(blicc_ld$polSs[blicc_ld$sp_i[sel_indx]:blicc_ld$sp_e[sel_indx]])))
    blicc_ld$polSs[blicc_ld$sp_i[sel_indx]:blicc_ld$sp_e[sel_indx]] <- 1.5  #default
  return(blicc_ld)
}



#' Sets gear selectivity mixture weights
#'
#' A single gear's mixture weights can be set or estimated if not provided. The
#' weights are given as positive values greater than zero. All weights are
#' relative to some base function which is set at 1.0, with other functions
#' having weights 0-1.
#'
#' If `mix_wt` is not provided it is estimated using a simple weighted least
#' squares procedure. If this is used, the selectivity parameters need to have
#' been previously set using [blip_set_sel].
#'
#' @export
#' @inheritParams blicc_mpd
#' @param gear Single integer or exact name indexing a gear
#' @param mix_wt Vector of prior mixture weights are required for each
#'   selectivity function after the first (which has a default weight of 1.0).
#'   If it is not provided it is loosely estimated from the data using weighted
#'   least squares.
#' @return The data object `blicc_ld` but with the new selectivity priors set
#'
blip_mix_wt <- function(blicc_ld,
                        gear,
                        mix_wt = NULL) {

  gear <- parse_gear(gear, blicc_ld)

  if (length(gear) != 1)
    stop("Error: `gear` must reference a single gear. \n")

  si <- 2L*gear-1L
  if (blicc_ld$GSmix1[si]==0)
    stop("Error: Gear ", as.character(gear), " has no mixtures. \n")

  mix_sel <-  blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]

  if (! is.null(mix_wt)) {
    if (! (is.vector(mix_wt) & is.numeric(mix_wt)) |
          length(mix_wt) != length(mix_sel) |
          any(mix_wt <= 0))
      stop("Error: `mix_wt` must be positive numeric values of length ", as.character(nmix), ". \n")
    blicc_ld$polSm[blicc_ld$NP + mix_sel] <- log(mix_wt)
  } else {
    sel_indx <- c(blicc_ld$GSbase[gear], blicc_ld$GSmix2[mix_sel])

    if (any(is.na(blicc_ld$polSm[blicc_ld$sp_i[sel_indx]])))
      stop("Error: Selectivity priors must be set using `blip_set_sel` before mixture weights can be estimated. \n")

    ns <- length(sel_indx)
    wts <- 1/(blicc_ld$fq[[gear]]+1) # approximate least-squares weights
    xnam <- paste0("x", sel_indx)
    model_df <- data.frame(y=pfq)
    model_df[, xnam] <- 0

    Sm <- exp(blicc_ld$polSm)
    for (i in 1:length(sel_indx)) {  # for each selectivity function in the mixture
      si <- sel_indx[i]
      switch(blicc_ld$fSel[sel_indx[i]],
             { #logistic- based very approx on integral of standard logistic F(x)=log(1+exp(x))
               model_df[[i+1L]] <- Rsel_logistic(Sm[blicc_ld$sp_i[si]:blicc_ld$sp_e[si]],  blicc_ld$LMP)
             },
             { #normal
               model_df[[i+1L]] <- Rsel_normal(Sm[blicc_ld$sp_i[si]:blicc_ld$sp_e[si]],  blicc_ld$LMP)
             },
             { #ssnormal
               model_df[[i+1L]] <- Rsel_ssnormal(Sm[blicc_ld$sp_i[si]:blicc_ld$sp_e[si]],  blicc_ld$LMP)
             },
             { #dsnormal
               model_df[[i+1L]] <- Rsel_dsnormal(Sm[blicc_ld$sp_i[si]:blicc_ld$sp_e[si]],  blicc_ld$LMP)
             }
      )

    }

    fmla <- stats::as.formula(paste("y ~ 0 + ", paste(xnam, collapse= "+")))
    ls <- stats::lm(fmla,  data=model_df, weights=wts)
    if (any(ls$coef < 0, na.rm = TRUE) | any(is.na(ls$coef))) {
      stop(paste0("Selectivity function for gear ",
                  blicc_ld$gear_name[gear],
                  " mixture weight estimation failed. Set the weight manually using `mix_wt`. \n"))
    }

    mix_wt <- ls$coef
    max_i <- which.max(mix_wt)
    mix_indx <- blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]
    if (max_i != 1L) { # switch main selectivity to GSbase
      blicc_ld$GSbase[gi] <- sel_indx[max_i]
      blicc_ld$GSmix2[mix_indx] <- sel_indx[-max_i]
    }
    blicc_ld$polSm[blicc_ld$NP+mix_indx] <- log(mix_wt[-max_i]/mix_wt[max_i])
  }

  return(blicc_ld)
}





#' Old: Sets data object from [blicc_dat()] to have a new selectivity priors
#'
#' Backup version of blip_selectivity.
#' The priors for each function are defined loosely based on the available data
#' ("empirical Bayes"). It is recommended that priors are weakly informative, so
#' they are set primarily to aid fitting and discourage values outside a
#' reasonable range. This estimate is based on trials for setting the
#' selectivity mode at each bin value and the steepness parameter (dependent on
#' the function) based on the 10% and 90% quartiles of the cumulative frequency
#' distribution. This provides a robust estimate of reasonable priors even with
#' selectivity mixtures. Alternatively, parameters can be set manually for each
#' gear using the sel_par and mix_wt arguments.
#'
#' For example, the location parameter is generally set with a log-normal mean
#' at the length frequency modes. Estimates far from this point would probably
#' have to be rejected during review after fitting anyway, so encouraging
#' estimates of the selectivity location parameters close to modes seems
#' reasonable. Prior parameters are selectivity function specific. Setting the
#' prior log-normal sd hyper-parameters polSs at values greater than the default
#' 1.5 will further reduce the influence of the prior setting. However, the
#' model will fit poorly if there is insufficient support in the data to fit the
#' parameters. In this case, the prior will determine the estimate.
#'
#' @inheritParams blicc_mpd
#' @inheritParams blicc_gear_sel
#' @param sel_par Vector of prior selectivity parameter means manually set for
#'   the model. Only one gear at a time can be specified to use this. Parameters
#'   must be in the correct order: in function order, with each function having
#'   a location parameter followed by one or two steepness parameters.
#' @param mix_wt Vector of mixture weights (0.0-1.0) are required for each
#'   selectivity function after the first (which has a default weight of 1.0).
##' @return The data object blicc_ld but with the new selectivity priors
#' @noRd
blip_selectivity_old <- function(blicc_ld,
                             gear = NULL,
                             sel_par = NULL,
                             mix_wt = NULL) {

  rec_mix_prior <- function(mix_i,      # current selectivity index
                            model_df,   # y and x values for model
                            mix_sel) {    # selectivity list

    # Recursive function to estimate prior parameters for mixtures as well
    # as single functions. Applies a grid search across bins for the location
    # parameter. This is an empirical Bayes approach.
    nx <- length(mix_sel)
    cur <- list(score=Inf, idx=integer(nx))
    ml <- mix_sel[[mix_i]]
    if (mix_i == nx) { # Calculate the least squares weights
      for (i in 1:ncol(ml)) {
        model_df[[mix_i+1L]] <- ml[,i]
        ls <- stats::lm(fmla,  data=model_df, weights=wts)
        ssq <- sum(ls$residuals^2)
        if (ssq < cur$score & ! any(ls$coef < 0, na.rm = TRUE) & ! any(is.na(ls$coef))) {
          cur$score <- ssq
          cur$idx[mix_i] <- i
          cur$wt <- ls$coef
        }
      }
    } else { # recurse
      for (i in 1:ncol(ml)) {
        model_df[[mix_i+1L]] <- ml[,i]
        test <- rec_mix_prior(mix_i+1L, model_df, mix_sel)
        if (test$score < cur$score) {
          cur <- test
          cur$idx[mix_i] <- i
        }
      }
    }
    return(cur)
  } # rec_mix_prior

  if (is.null(gear))
    gear <- 1:blicc_ld$NG
  else
    gear <- parse_gear(gear, blicc_ld)


  if (! is.null(sel_par)) {
    if (length(gear > 1)) stop("Error: A single gear must be specified if sel_par is used. \n")
    si <- 2L * gear - 1L
    sel_indx <- blicc_ld$GSbase[gear] # base function index

    if (blicc_ld$GSmix1[si] > 1) {
      nm <- with(blicc_ld, GSmix1[si+1L] - GSmix1[si] + 1L) # number of mixtures
      if (nm != length(mix_wt) | any(mix_wt <= 0))
        stop(paste0("Error: ", nm,
                    " prior mixture weights need to be specified for this gear. \n"))
      mix_indx <- with(blicc_ld, GSmix2[GSmix1[si]:GSmix1[si+1L]])

      blicc_ld$polSm[blicc_ld$NP+mix_indx] <-  log(mix_wt) # default
      blicc_ld$polSs[blicc_ld$polSs==0 & c(rep(FALSE, blicc_ld$NP), rep(TRUE, blicc_ld$NM))] <- 1.5 # default
      sel_indx <- c(sel_indx, mix_indx)  # indices for mixture functions
    }

    par_indx <- integer(0)
    for (i in sel_indx) {
      par_indx <- with(blicc_ld, c(par_indx, sp_i[i]:sp_e[i]))
    }
    if (length(par_indx) != length(sel_par) | any(sel_par <= 0))
      stop(paste0("Error: length of sel_par must be positive and ",
                  "equal number of parameters for the gear (",
                  length(par_indx), "). \n"))
    blicc_ld$polSm[par_indx] <- log(sel_par)
    blicc_ld$polSs[par_indx] <- 1.5 # default
    return(blicc_ld)
  } # if (! is.null(sel_par))

  # else

  ssd <- 1.281552
  ssd <- 2.0
  Galpha <- exp(blicc_ld$polGam)
  Linf <- blicc_ld$poLinfm
  Gbeta <- Galpha/Linf
  LLB <- blicc_ld$LLB
  LMP <- blicc_ld$LMP
  Zk <- exp(blicc_ld$polMkm)*blicc_ld$M_L
  gl <- statmod::gauss.quad(110, "laguerre", alpha=0)
  pop <- Rpop_len(gl$nodes, gl$weights, LLB, Zk, Galpha, Gbeta)

  for (gi in gear) {
    pfq <- blicc_ld$fq[[gi]]/pop  # adjust data for mortality
    pfq <- pfq/sum(pfq)           # normalise
    cfq <- cumsum(pfq)            # cumulative sum
    i01 <- max(which(cfq<0.1))    # 10% quartile index
    i05 <- min(which(cfq>0.5))    # 50% quartile index
    i09 <- min(which(cfq>0.9))    # 90% quartile index
    si <- 2L * gi - 1L
    if (blicc_ld$GSmix1[si] > 0L)
      ns <- blicc_ld$GSmix1[si+1L] - blicc_ld$GSmix1[si] + 2L
    else
      ns <- 1L

    start_index <- 3L
    end_index <- blicc_ld$NB - start_index
    win_len <- (end_index - start_index) %/% ns
    if (win_len <= 6) stop(paste0("Error: Too many mixture functions for the number of bins: gear ",
                                  blicc_ld$gear_name[gi]))
    lowidth <- win_len %/% 2
    hiwidth <- win_len - lowidth

    sel_indx <- blicc_ld$GSbase[gi]
    if (ns > 1) sel_indx <- with(blicc_ld, c(sel_indx, GSmix2[GSmix1[si]:GSmix1[si+1L]]))
    sel_fun <- blicc_ld$fSel[sel_indx]

    if (ns > 3) {
      # in this case equally spaced windows each with a cumulative calculation to
      # set parameters
      # issue warning
      warning("More than 3 mixture functions is not recommended.
              Check the priors. You may need to set the selectivity parameter
              priors manually. \n")
      LMP <- blicc_ld$LMP
      index <- rep(1:ns, each=win_len)
      if (length(index) < (blicc_ld$NB-6L)) index <- c(index, rep(ns, (blicc_ld$NB-6L)-length(index)))
      index <- c(rep(1L, start_index), index, rep(ns, start_index))
      sum_pfq <- unlist(tapply(pfq, index, sum))
      cfq <- tapply(cfq, index, FUN=function(x) x/max(x))
      ssd <- 0.6745
      base_index <- cumsum(c(0, sapply(cfq, length)[-ns]))
      lo <- LMP[sapply(cfq, function(x) max(1, which(x<0.25)), USE.NAMES = FALSE) + base_index] - 0.5
      loc <- LMP[sapply(cfq, function(x) min(which(x>0.5)), USE.NAMES = FALSE) + base_index]
      hi <- LMP[sapply(cfq, function(x) min(which(x>0.75), length(x)), USE.NAMES = FALSE) + base_index] + 0.5
      # move the most significant frequency window to the base
      max2min <- order(sum_pfq, decreasing=T)
      lo <- lo[max2min]
      loc <- loc[max2min]
      hi <- hi[max2min]
      sum_pfq <- sum_pfq[max2min]  # empirical basis for the mixture weights
      for (i in 1:ns) {
        switch(sel_fun[i],
               { #logistic
                 blicc_ld$polSm[blicc_ld$sp_i[sel_indx[i]]:blicc_ld$sp_e[sel_indx[i]]] <-
                   log(c(loc[i], -0.1*log(1+exp(0.2*(hi[i]-loc[i])))/(lo[i]-loc[i])))
               },
               { #normal
                 blicc_ld$polSm[blicc_ld$sp_i[sel_indx[i]]:blicc_ld$sp_e[sel_indx[i]]] <-
                   log(c(loc[i], (0.5*(hi[i]-lo[i])/ssd)^-2))
               },
               { #ssnormal
                 blicc_ld$polSm[blicc_ld$sp_i[sel_indx[i]]:blicc_ld$sp_e[sel_indx[i]]] <-
                   log(c(loc[i], ((lo[i]-loc[i])/ssd)^-2))
               },
               { #dsnormal
                 blicc_ld$polSm[blicc_ld$sp_i[sel_indx[i]]:blicc_ld$sp_e[sel_indx[i]]] <-
                   log(c(loc[i], ((loc[i]-lo[i])/ssd)^-2, ((hi[i]-loc[i])/ssd)^-2))
               }
        )
      }
      mix_indx <- blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]
      blicc_ld$polSm[blicc_ld$NP+mix_indx] <- log((sum_pfq[-1]+1)/(sum_pfq[1]+1))
    } else {
      # 1-3 selectivity functions, so we try to "fit" them

      # calculate spread for each location using windows at 25% and 75%iles.
      #   Set lower/upper boundary at effective zero if outside the window
      # calculate selectivity using each location and spread for each mixture function

      sel_mix <- list()
      sel_par_list <- list()
      base_loc <- LMP[start_index:end_index] # range of location parameters
      lo <- numeric(length(base_loc))
      hi <- lo
      for (j in 1:length(base_loc)) {
        win <- (start_index+j-lowidth):(start_index+j+hiwidth)
        win <- win[win >= 1 & win <= length(LMP)]
        cs <- cumsum(cfq[win])
        cs <- cs/cs[length(cs)]
        lo[j] <- LMP[max(1, which(cs<0.25))] - 0.5
        hi[j] <- LMP[min(which(cs>0.75))] + 0.5
      }

      for (i in 1:ns) {  # for each selectivity function in the mixture
        switch(sel_fun[i],
               { #logistic- based very approx on integral of standard logistic F(x)=log(1+exp(x))
                 spar <- cbind(base_loc,
                               abs(0.5*(log(exp(0.1)-1)/(lo-base_loc) + log(exp(0.9)-1)/(hi-base_loc))))
                 sel_mat <- apply(spar, 1, Rsel_logistic, LMP=LMP)
               },
               { #normal
                 if (i==1) spar <- cbind(base_loc, ((lo-base_loc)/ssd)^-2)
                 else if (i==ns) spar <- cbind(base_loc, ((hi-base_loc)/ssd)^-2)
                 else spar <- cbind(base_loc, (0.5*(hi-lo)/ssd)^-2)
                 sel_mat <- apply(spar, 1, Rsel_normal, LMP=LMP)
               },
               { #ssnormal
                 spar <- cbind(base_loc, ((lo-base_loc)/ssd)^-2)
                 sel_mat <- apply(spar, 1, Rsel_ssnormal, LMP=LMP)
               },
               { #dsnormal
                 spar <- cbind(base_loc, ((base_loc-lo)/ssd)^-2, ((hi-base_loc)/ssd)^-2)
                 sel_mat <- apply(spar, 1, Rsel_dsnormal, LMP=LMP)
               }
        )
        sel_mix[[i]] <- sel_mat
        sel_par_list[[i]] <- spar
      }

      # Set up fit also used for mixture weights
      wts <- 1/(blicc_ld$fq[[gi]]+1) # approximate least-squares weights
      xnam <- paste0("x", 1:ns)
      model_df <- data.frame(y=pfq)
      model_df[, xnam] <- 0
      fmla <- stats::as.formula(paste("y ~ 0 + ", paste(xnam, collapse= "+")))

      min_res <- rec_mix_prior(1L, model_df, sel_mix)

      if (any(min_res$idx==0)) stop(paste0("Selectivity function for gear ",
                                           blicc_ld$gear_name[gi],
                                           " prior estimation failed. \n"))

      for (i in 1L:ns) {
        blicc_ld$polSm[blicc_ld$sp_i[sel_indx[i]]:blicc_ld$sp_e[sel_indx[i]]] <-
          log(sel_par_list[[i]][min_res$idx[i],])
      }

      if (ns > 1L) {
        max_i <- which.max(min_res$wt)
        mix_indx <- blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]
        if (max_i != 1L) { # switch main selectivity to GSbase
          blicc_ld$GSbase[gi] <- sel_indx[max_i]
          blicc_ld$GSmix2[mix_indx] <- sel_indx[-max_i]
        }
        blicc_ld$polSm[blicc_ld$NP+mix_indx] <- log(min_res$wt[-max_i]/min_res$wt[max_i])
      }
    }
  }
  blicc_ld$polSs[] <- 1.5  # default

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
