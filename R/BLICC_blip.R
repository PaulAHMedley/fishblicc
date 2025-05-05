# BLICC Setting Priors functions ----------------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><



#' Set `Linf` prior parameters in the data object
#'
#' New `Linf` prior hyper-parameters are set in the [blicc_dat] data object as a
#' vector of the mean (mu) and standard deviation (sigma) for the normal prior
#' to be used. The `Linf` will need to be a vector of 2 values.
#'
#' The `Linf` defines the gamma growth probability density function mean, so it
#' is the mean asymptotic length of the population.
#'
#' @export
#' @inheritParams blicc_mpd
#' @param Linf  A numeric vector of double containing the mu and sigma for the
#'   prior normal.
#' @param model_name A string for a replacement model name in the data object.
#'   Optional.
#' @return The data object blicc_ld but with the prior for `Linf` changed.
#' @examples
#' new_ld <- blip_Linf(gillnet_ld, c(60,2), model_name="Sensitivity")
#' 
blip_Linf <- function(blicc_ld,
                      Linf,
                      model_name=NULL) {
  if (! is.vector(Linf, mode="double") | length(Linf) != 2)
    stop("Error: Linf must be a vector of 2 (mu, sigma) for the prior. \n")
  if (Linf[1] <= min(blicc_ld$LLB)) {
    stop(paste(
      "Error: Linf must be greater than the lowest bin value:",
      as.character(min(blicc_ld$LLB)), " \n"
    ))
  }
  if (Linf[2] <= 0) {
    stop("Error: Linf prior sd must be greater than zero. \n")
  }

  if (!is.null(model_name))
    blicc_ld$model_name <- model_name
  # Expected Linf
  blicc_ld$poLinfm <- Linf[1]
  # sd for the normal Linf, see above
  blicc_ld$poLinfs <- Linf[2]
  return(blicc_ld)
}


#' Set `Galpha` log-normal prior hyper-parameters in the data object
#'
#' The `Galpha` lognormal prior is updated with new values for the mean (mu) and
#' standard deviation (sigma). The CV is likely to be between 5% (`Galpha` mu =
#' `log(1/0.05^2)`) and 30% (`Galpha` mu = `log(1/0.3^2)`) Values outside this
#' range are not recommended. A 30% CV makes length very uninformative on age.
#'
#' `Galpha` is the alpha parameter in the gamma distribution for the Linf
#' probability density. The default in [blicc_dat] is CV=10%: `lGa1pha = c(log(1
#' / 0.1 ^ 2), 0.25)`. 
#'
#' @export
#' @inheritParams blip_Linf
#' @param lGalpha  A numeric vector of double containing the mean and sd for the
#'   prior log-normal.
#' @return The data object blicc_ld but with the prior for Galpha changed.
#' @examples
#' new_ld <- blip_Galpha(gillnet_ld, lGalpha=c(log(1/0.05^2), 0.1))
#' 
blip_Galpha <- function(blicc_ld,
                        lGalpha) {
  if (!(is.numeric(lGalpha) & length(lGalpha)==2 & lGalpha[2] > 0))
    stop("Error: LGalpha must be numeric vector of the mu,sd for the Galpha lognormal prior. \n")
  if (lGalpha[1] > log(1/0.05^2) | lGalpha[1] < log(1/0.3^2))
    warning("The Galpha prior lognormal mean is outside the expected 5-30% CV. \n")
  # Growth mean CV
  blicc_ld$polGam <- lGalpha[1]
  blicc_ld$polGas <- lGalpha[2]
  return(blicc_ld)
}


#' Set `Mk` log-normal prior hyper-parameters in the data object
#'
#' The `Mk` hyper-parameters are provided as a mean and sigma for the lognormal.
#' The function checks the proposed values are valid. Also, if the ref_length
#' parameter is defined, apply the inverse length function for natural
#' mortality. These are replaced in the data object, which is then returned.
#' Note that the Mk must be provided as the log value. The default in
#' [blicc_dat] depends on the length at 50% maturity, but might be reasonably
#' close to: `lMk = c(log(1.5), 0.1)`. Values are only changed if the supplied
#' value is not `NA`.
#'
#' The `ref_length` parameter, if used, indicates a length-inverse natural
#' mortality model, where `ref_length` is the length where the given natural
#' mortality applies. Therefore, the natural mortality within each length bin
#' will be calculated as `Mki = Mk*ref_length/L` where `L` is the mid-point for
#' the i^th^ bin. Natural mortality will be constant within each bin. If
#' `ref_length` is not given, a constant natural mortality is assumed.
#'
#' @export
#' @inheritParams blip_Linf
#' @inheritParams blicc_dat
#' @param lMk A vector of the mu and sigma for the lognormal natural mortality
#'   prior
#' @param ref_length Reference length in the length-inverse mortality is
#'   applied. Set to -1 to turn off the length-inverse model.
#' @return The supplied data object blicc_ld but with the prior and function for
#'   Mk changed.
#' @examples
#' new_ld <- blip_Mk(gillnet_ld, lMk=c(log(1.7), NA), ref_length=25)
#' 
blip_Mk <- function(blicc_ld,
                    lMk = c(NA_real_, NA_real_),
                    ref_length = -1,
                    model_name = NULL) {
  # Natural mortality
  if (ref_length > 0){
    if (ref_length < min(blicc_ld$LLB) | ref_length > max(blicc_ld$LLB))
      stop("Error: The reference length for the natural mortality must be within the length frequencies. \n")
    M_L <- ref_length/blicc_ld$LMP
  } else {
    M_L <- rep(1, blicc_ld$NB)  # Fixed natural mortality
    ref_length <- -1
  }
  blicc_ld$M_L <- M_L
  blicc_ld$ref_length <- ref_length

  if (! (is.vector(lMk, mode = "numeric") & length(lMk)==2))
    stop("Error: natural mortality must be provided as vector of mean and sigma for the lognormal. \n")

  if (!is.null(model_name))
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


#' Set `Fk` log-normal prior hyper-parameters in the data object
#'
#' The `Fk` lognormal prior is updated with new mu for each gear and a single
#' sigma parameter. If no value (`NA`) is provided for `lFk`, the defaults are
#' assigned. These are the current natural mortality for fishing mortality mu
#' and 2.0 for the sigma parameter. The sigma hyper-parameter should be kept
#' reasonably large to allow the model to estimate this parameter. In general,
#' the default values should be sufficient as the prior should be
#' non-informative.
#'
#' For the multigear model, if a gear has negligible catches so the relative
#' catch vector indicates zero catch for that gear, no fishing mortality is
#' estimated for it, and the gear's `Fk` prior is skipped as it is fixed at
#' zero.
#'
#' @export
#' @inheritParams blip_Linf
#' @param lFk   A vector containing the log-normal mean `Fk` for each gear
#'   having a non-negligible catch
#' @param lFks  A double containing the same log-normal sigma `Fk` for all gears
#' @return The data object `blicc_ld` with the prior for `Fk` changed.
#' @examples
#' new_ld <- blip_Fk(gillnet_ld, lFk=log(1.9), lFks=1.5)
#' 
blip_Fk <- function(blicc_ld,
                    lFk  = NULL,
                    lFks = NULL,
                    model_name=NULL) {
  if (is.null(lFk)) {
    lFk <- blicc_ld$polMkm + log(blicc_ld$prop_catch)
  }
  if (any(is.null(lFks)))
    lFks = 2.0
  if (length(lFks) != 1) 
    stop("Error: lFks must be a single real number for the prior Fk log-normal sigma.")
  if (! (is.numeric(lFk) & is.numeric(lFks)))
    stop("Error in blip_Fk: supplied values are not numeric.")
  if (blicc_ld$NF != length(lFk))
    stop(paste("Error in blip_Fk: lFk array not of length: ", as.character(blicc_ld$NF)))
  if (! is.null(model_name))
    blicc_ld$model_name <- model_name
  blicc_ld$polFkm <- as.array(lFk)
  blicc_ld$polFks <- lFks
  return(blicc_ld)
}


#' Set life-history parameters and vectors in the data object
#'
#' Sets the maturity-at-length and weight-at-length vectors in the data object
#' (from [blicc_dat]) to new values. The values are either provided as vectors
#' or model parameters to calculate the vectors. Calculations are done for the
#' length bin mid-points. The models for the calculation are the standard
#' length-weight (`W=aL^b`) and logistic model for the maturity-at-length: (`Lm
#' = 1/(1+exp(Ls*(L-L50)))`).
#'
#' @export
#' @inheritParams blip_Linf
#' @inheritParams blicc_dat
#' @param set_defaults Logical indicating whether to set defaults if parameters
#'   are `NA`. Leaves them alone if `FALSE`.
#' @return The data object blicc_ld with the new life-history parameters
#' @examples
#' new_ld <- blip_LH(gillnet_ld, a=1.2e-4, b=3, model_name="Alternative LW")
#' 
blip_LH <-
  function(blicc_ld,
           a = NULL,
           b = NULL,
           L50 = NULL,
           L95 = NULL,
           ma_L = NULL,
           wt_L = NULL,
           model_name = NULL,
           set_defaults = FALSE) {
    Linf <- blicc_ld$poLinfm
    # Weight and maturity
    if (is.null(L50)) {
      if (set_defaults)
        L50 <- 0.66 * Linf
    } else {
      if (L50 <= 0.2 * Linf | L50 >= Linf) {
        stop("Error: Length at 50% maturity must be greater than 0.2*Linf and less than Linf. \n")
      }
    }
    if (is.null(L95)) {
      # Ls = -log(1/0.95 - 1)/(L95-L50)
      if (set_defaults)
        Ls <- -log(1 / 0.95 - 1) / (0.05 * (Linf - L50))
      else
        Ls <- NULL
    } else {
      if (L95 <= L50) {
        stop("Error: Length at 95% maturity must be greater than L50. \n")
      }
      if (L95 > Linf) {
        warning("Length at 95% maturity is greater than Linf. \n")
      }
      Ls <- -log(1 / 0.95 - 1) / (L95 - L50)
    }
    
    if (is.null(wt_L)) {
      if (is.null(a) & set_defaults) {
        warning("No weight-at-length information provided - the weight units will be incorrect. \n")
        a <- 1.0
      }
      if (is.null(b)) {
        if (set_defaults) b <- 3.0 
      } else if (! is.na(b)) {
        if (b <= 2 | b > 4) 
          stop("Error: Length-weight exponent (b) must be greater than 2 and less than 4. \n")
      }
      if (!(is.null(a) | is.null(b)))
        wt_L <-
          with(blicc_ld, a * exp(b * log(LMP)))    # Estimated biomass per recruit
      else {
        if (set_defaults)
          warning("a or b not specified: weight-at-length not changed. \n")
      }
    } else {
      if (!is.numeric(wt_L) | length(wt_L) != blicc_ld$NB) {
        stop("Error: Weight-at-length vector must be a numeric vector with size equal to the number of length bins. \n")
      }
    }
    
    if (any(is.null(ma_L))) {
      if (!(is.null(Ls) | is.null(L50)))
        ma_L <-
          with(blicc_ld, wt_L / (1 + exp(-Ls * (LMP - L50))))    #Mature biomass
      else
        warning("L50 or L95 not specified: mature biomass -at-length not changed. \n")
    } else {
      if (!is.numeric(ma_L) | length(ma_L) != blicc_ld$NB) {
        stop(
          "Error: Mature biomass -at-length vector must be a numeric vector with size equal to the number of length bins. \n"
        )
      }
    }
    
    if (!is.null(model_name))
      blicc_ld$model_name <- model_name
    if (!is.null(a))
      blicc_ld$a <- a
    if (!is.null(b))
      blicc_ld$b <- b
    if (!is.null(L50))
      blicc_ld$L50 <- L50
    if (!is.null(Ls))
      blicc_ld$Ls <- Ls
    if (!any(is.null(wt_L)))
      blicc_ld$wt_L <- wt_L
    if (!any(is.null(ma_L)))
      blicc_ld$ma_L <- ma_L
    return(blicc_ld)
  }

#'Set automatic selectivity prior parameters in data object
#'
#'The priors for each function are defined loosely based on the available data
#'("empirical Bayes"). It is recommended that priors are weakly informative, so
#'they are set primarily to aid fitting and discourage values outside a
#'reasonable range.
#'
#'The prior hyper-parameters are "estimated" if there is a single selectivity
#'for a gear. The estimates are based on the 50% (median), 10% and 90% quartiles
#'of the cumulative frequency. This ensures that the selectivity function is
#'centred on the frequency data and provides a robust estimate of reasonable
#'priors. Alternatively, setting parameters manually can be done using the
#'[blip_sel] function.
#'
#'If a selectivity is made up of a mixture of functions, prior parameters are
#'not estimated and [blip_sel] must be used to set the parameters manually. This
#'is necessary because the user must propose the hypothesis for the mixtures and
#'it cannot easily be determined from the data what a suitable hypothesis might
#'be.
#'
#'@export
#'@inheritParams blicc_selfun
#'@return The data object blicc_ld but with the new selectivity priors
#' @examples
#'new_ld <- blip_sel_auto(gillnet_ld)
#'
blip_sel_auto <- function(blicc_ld,
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
  mort_corr <- 1.0/(pop + (1.0/blicc_ld$NB)) # robust mortality correction
  
  for (si in sel_indx) {
    par_range <- with(blicc_ld, sp_i[si]:sp_e[si])
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
          "See function `blip_sel`. \n"
        )
      )
    else {
      
      qi <- which(blicc_ld$Gi==gi) # combine frequencies for this gear
      fq <- Reduce(`+`, blicc_ld$fq[qi])
      pfq <- fq * mort_corr  # adjust data for mortality
      pfq <- pfq / sum(pfq)                 # normalise
      cfq <- cumsum(pfq)                    # cumulative sum
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


#'Set a selectivity function's prior hyper-parameters to the given values
#'
#'Each parameter set is for a single selectivity function indexed by `sel_indx`.
#'The location parameter (`loc`) is either the 50% selectivity for the logistic
#'or the mode of the particular normal selectivity function, and is required.
#'
#'The `lslope` parameter is a single or double log value of the slope for the
#'relevant function. This is either the log of the logistic slope parameter, or
#'the log of the reciprocal of the variance for the normal functions. For the
#'double-sided normal, two parameters must be provided and for all other
#'functions only one parameter. If `lslope` is not provided, a default slope of
#' -4.5 is applied if the slope parameters have not already been set. If values are
#'already present these are conserved. Mixture weights, if they are used, can be
#'set in [blip_mix].
#'
#'@export
#'@inheritParams blicc_mpd
#'@param sel_indx Single integer indexing a selectivity function
#'@param loc Single location parameter for the selectivity function
#'@param lslope Vector of prior selectivity log slope parameters manually set
#'  for the specified selectivity function. Optional.
#'@return The data object `blicc_ld` but with the new selectivity priors set
#'@examples
#'new_ld <- blip_sel(gillnet_ld, sel_indx = 1, loc = 30, lslope = -2.8)
#'
blip_sel <- function(blicc_ld,
                     sel_indx,
                     loc,
                     lslope = NULL) {
  
  sel_indx <- parse_sel_indx(sel_indx, blicc_ld, TRUE)

  if (! (is.vector(loc, mode="numeric")) |
      length(loc) != 1 |
      loc <= 0)
    stop("Error: `loc` must be a single positive numeric value for logistic 50% selectivity or normal mode. \n")

  if (blicc_ld$fSel[sel_indx] == 4) np <- 2 else np <- 1

  if (! is.null(lslope)) {
    if (! (is.vector(lslope, mode="numeric")) |
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



#' Set gear selectivity mixture weights
#'
#' A single gear's prior mixture weights can be set or estimated if not
#' provided. The weights are given as positive values greater than zero. All
#' weights are relative to some base function which is set at 1.0, with other
#' functions having weights 0-1.
#'
#' If `mix_wt` is not provided it is estimated using a simple weighted least
#' squares procedure. If this is used, the selectivity parameters for each
#' function will need to have been previously set using [blip_sel].
#'
#' @export
#' @inheritParams blicc_mpd
#' @param gear Single integer or exact name indexing one gear
#' @param mix_wt Vector of prior mixture weights are required for each
#'   selectivity function after the first (which has a default weight of 1.0).
#'   If it is not provided it is loosely estimated from the data using weighted
#'   least squares.
#' @return The data object `blicc_ld` but with the new selectivity priors set
#' 
blip_mix <- function(blicc_ld,
                     gear,
                     mix_wt = NULL) {
  nmix=pfq=gi=NULL
  
  gear <- parse_gear(gear, blicc_ld)

  if (length(gear) != 1)
    stop("Error: `gear` must reference a single gear. \n")

  si <- 2L*gear-1L
  if (blicc_ld$GSmix1[si]==0)
    stop("Error: Gear ", as.character(gear), " has no mixtures. \n")

  mix_sel <-  blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]

  if (! is.null(mix_wt)) {
    if (! (is.vector(mix_wt, mode="numeric")) |
          length(mix_wt) != length(mix_sel) |
          any(mix_wt <= 0))
      stop("Error: `mix_wt` must be positive numeric values of length ", as.character(nmix), ". \n")
    blicc_ld$polSm[blicc_ld$NP + mix_sel] <- log(mix_wt)
  } else {
    sel_indx <- c(blicc_ld$GSbase[gear], blicc_ld$GSmix2[mix_sel])

    if (any(is.na(blicc_ld$polSm[blicc_ld$sp_i[sel_indx]])))
      stop("Error: Selectivity priors must be set using `blip_sel` before mixture weights can be estimated. \n")

    ns <- length(sel_indx)
    wts <- 1/(blicc_ld$fq[[gear]]+1) # approximate least-squares weights
    xnam <- paste0("x", sel_indx)
    model_df <- data.frame(y=pfq)
    model_df[, xnam] <- 0

    Sm <- exp(blicc_ld$polSm)
    for (i in 1:length(sel_indx)) {  # for each selectivity function in the mixture
      si <- sel_indx[i]
      switch(blicc_ld$fSel[si],
             { 
               model_df[[i+1L]] <- with(blicc_ld, Rsel_logistic(Sm[sp_i[si]:sp_e[si]],  LMP))
             },
             { #normal
               model_df[[i+1L]] <- with(blicc_ld, Rsel_normal(Sm[sp_i[si]:sp_e[si]], LMP))
             },
             { #ssnormal
               model_df[[i+1L]] <- with(blicc_ld, Rsel_ssnormal(Sm[sp_i[si]:sp_e[si]], LMP))
             },
             { #dsnormal
               model_df[[i+1L]] <- with(blicc_ld, Rsel_dsnormal(Sm[sp_i[si]:sp_e[si]],  LMP))
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


#' Normalize selectivity mixture weights, so the maximum mixture weight is 1.0
#'
#' The function checks each gear that has mixtures that the maximum weight is
#' 1.0. If it is higher than 1.0, the selectivity with the highest weight is
#' moved to the base component using `blip_MainComp`.
#' 
#' @export
#' @inheritParams blicc_mpd
#' @return The data object `blicc_ld` with normalized weights
#' 
blip_normalize_mix <- function(blicc_ld) {
  for (gi in 1:blicc_ld$NG) {
    si <- 2L*gi-1L
    if (blicc_ld$GSmix1[si] > 0) {
      mix_indx <-  blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]
      max_i <- which.max(blicc_ld$polSm[blicc_ld$NP+mix_indx])
      if (blicc_ld$polSm[blicc_ld$NP+mix_indx[max_i]] > 0) {
        blicc_ld <- blip_main_sel(blicc_ld, gi, blicc_ld$GSmix2[mix_indx[max_i]])
      }
    }
  }
  return(blicc_ld)
}


#' Set gear selectivity component to main
#'
#' A gear's selectivity component is set so that its weight equals 1.0. Other 
#' mixture weights are adjusted accordingly. 
#'
#' @export
#' @inheritParams blip_mix
#' @param sel_N Index for the new main selectivity component with weight=1.
#' @return The data object `blicc_ld` but with the new main selectivity
#' 
blip_main_sel <- function(blicc_ld,
                          gear,
                          sel_N) {

  if (length(gear) != 1)
    stop("Error: `gear` must reference a single gear. \n")
  
  gear <- parse_gear(gear, blicc_ld)
  
  if (length(sel_N) != 1 | ! is.numeric(sel_N))
    stop("Error: `sel_N` must reference a single selectivity component. \n")
  
  # check component is in gear
  
  si <- 2L*gear-1L
  if (blicc_ld$GSmix1[si]==0) return(blicc_ld) # ignore as no mixtures

  mix_sel_i <-  blicc_ld$GSmix1[si]:blicc_ld$GSmix1[si+1L]
  new_sel_i <- blicc_ld$GSmix2[mix_sel_i]==sel_N
  
  if (! any(new_sel_i)) 
    stop("Error: Selectivity component not used by this gear. \n")

  old_sel <- blicc_ld$GSbase[gear]
  blicc_ld$GSbase[gear] <- sel_N
  sel_i <- mix_sel_i[new_sel_i]
  blicc_ld$GSmix2[sel_i] <- old_sel
  wt_norm <- blicc_ld$polSm[blicc_ld$NP+sel_i]
  blicc_ld$polSm[blicc_ld$NP+sel_i] <- 0
  blicc_ld$polSm[blicc_ld$NP+mix_sel_i] <- blicc_ld$polSm[blicc_ld$NP+mix_sel_i] - wt_norm

  return(blicc_ld)
}


#' Set the `NB_phi` prior hyper-parameters in the data object
#'
#' The `NB_phi` log-normal prior is updated with a mean (mu) and standard
#' deviation (sigma). NB_phi controls the negative binomial over-dispersion
#' compared to the Poisson distribution. It is sometimes useful to reduce the
#' NB_phi lognormal standard deviation hyper-parameter below the 0.5 default to
#' avoid variation in the data being interpreted as observation error when it
#' can be explained better by the model.
#'
#' The variance for the negative binomial is ` Var = mu + mu^2 / NB_phi`
#' compared to `Var = mu` for the Poisson distribution.
#'
#' @export
#' @inheritParams blip_Linf
#' @param lNBphi  A numeric vector of double containing the mu and sigma for the
#'   prior log-normal.
#' @return The data object blicc_ld but with the prior for `NBphi` changed.
#' @examples
#' new_ld <- blip_NBphi(gillnet_ld, lNBphi = c(log(100), 0.05))
#' 
blip_NBphi <- function(blicc_ld,
                       lNBphi) {
  if (!is.numeric(lNBphi) & length(lNBphi)==2)
    stop("Error: lNBphi must be numeric vector of the c(mu,sd) for the NB_phi lognormal prior. \n")
  blicc_ld$polNB_phim <- lNBphi[1]
  blicc_ld$polNB_phis <- lNBphi[2]
  return(blicc_ld)
}


#' Set the relative catch log-normal likelihood standard deviation
#'
#' The relative catches for multigear fisheries are fitted using a lognormal. The
#' default is set at around 5%, but lower values may be necessary to ensure the
#' catches are fitted particularly when large number of length observations are
#' used. Has no effect in single gear fisheries.
#'
#' @export
#' @inheritParams blip_Linf
#' @param sigma  The standard deviation for the lognormal likelihood
#' @return The data object blicc_ld but with the value for `polCs` changed.
#' @examples
#' new_ld <- blip_catchLN(trgl_ld, sigma = 0.01)
#'
blip_catchLN <- function(blicc_ld,
                         sigma) {
  if (!is.numeric(sigma) | length(sigma)>1 | any(sigma <= 0))
    stop("Error: sigma must be a single positive numeric value. \n")
  blicc_ld$polCs <- sigma
  return(blicc_ld)
}


#' Estimate a minimum number of Gauss-Laguerre quadrature nodes for a defined
#' tolerance
#'
#' The procedure compares population estimates at length to estimates using a
#' safe number of nodes to identify the smallest number of nodes that produce
#' estimates that are still within the specified tolerance.
#'
#' @inheritParams blicc_mpd
#' @param draws A dataframe of parameter draws to test
#' @param toler The minimum acceptable tolerance for the integral
#' @return A list of Gauss-Laguerre nodes and weights meeting the minimum
#'   tolerance
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
    for (qi in 1:blicc_ld$NQ) {
      if (blicc_ld$Fkq[qi] > 0) {
        Zki[[i]] <- Zki[[i]] + FSel[[blicc_ld$Gi[qi]]] * Fk[i, blicc_ld$Fkq[qi]] # Total mortality
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
      if (is.na(err)) stop("Error: Numerical error in LG integration: check the data. \n")
      else if (err > toler) break
    }
    if (err < toler | Opt_NK >= 110) break
    Opt_NK <- Opt_NK+10
  }
  return(Opt_NK)
}
