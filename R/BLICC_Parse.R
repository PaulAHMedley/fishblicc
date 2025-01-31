# BLICC Parsing Functions for internal use ------------------------------------

# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><
# ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <>< ><> <><


#' Return a list of supported selectivity functions with parameters
#'
#' This function documents the list of available selectivity functions. This is
#' used in various routines as a constant. The intent is to make it easier to
#' expand the number of functions in future. Note however that the method now
#' offers mixtures of selectivity functions which is very flexible, so this list
#' may not expand much.
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
#' The selectivity function parameter is converted to an integer index of the
#' function. The selectivity function can be referenced by its name or an
#' integer index.
#'
#' @inheritParams blicc_selfun
#' @return Integer gear index or NA if no gear can be specified
#' @noRd
#' 
parse_selectivity <- function(sel_fun, blicc_ld) {
  func <- Rsel_functions()
  Nfunc <- length(func$short_name)
  func_list <- paste(func$short_name, collapse=" ")
  errmsg <- paste0("Error: specified selectivity function must be the function name (",
                   func_list,
                   " or all) or an integer between 1 and ", as.character(Nfunc))
  if (is.character(sel_fun)){
    if (! all(sel_fun %in% func$short_name))
      stop(errmsg)
    sel_fun <- match(sel_fun, func$short_name)
  } else {
    if (! is.numeric(sel_fun))
      stop(errmsg)
    sel_fun <- as.integer(sel_fun)
    if (min(sel_fun) < 1 | max(sel_fun) > Nfunc)
      stop(errmsg)
  }
  return(sel_fun)
}


#' Parse the selectivity index, converting to an integer if necessary
#'
#' The selectivity index parameter is checked to be within bounds and is
#' converted to an integer index if necessary.
#'
#' @inheritParams selfun
#' @param SingleValue True or False, enforces a single value if required.
#' @return Integer index of the selectivity functions in the model (or stops
#'   with error message)
#' @noRd
#' 
parse_sel_indx  <- function(sel_indx, blicc_ld, SingleValue = FALSE) {
  if (SingleValue & (length(sel_indx) != 1))
    stop(paste0("Error: sel_indx must be an integer with value, when rounded, between 1 and ",
                     as.character(blicc_ld$NS), ". \n"))
  
  if (! (is.vector(sel_indx) & is.numeric(sel_indx)) |
      any(round(sel_indx) < 1) |
      any(round(sel_indx) > blicc_ld$NS))
    stop(paste0("Error: sel_indx must be integers which, when rounded, are between 1 and ",
                as.character(blicc_ld$NS), ". \n"))
  
  sel_indx <- unique(round(sel_indx))
  
  return(sel_indx)
}



#' Parse the gear parameter, converting to an integer
#'
#' The gear parameter is converted to an integer index of the gear. Exact gear
#' names are also accepted.
#'
#' @inheritParams blicc_mpd
#' @param Gear  A vector of 1 or more numbers or strings representing valid
#'   gears in the model
#' @return Integer gear index (or stops with error message)
#' @noRd
#' 
parse_gear <- function(Gear, blicc_ld) {
  if (blicc_ld$NG == 1) {
    return(1)
  } else {
    if (any(is.na(Gear))) {
      stop(
        paste0(
          "Error: Gears must be specified as 'All', or exact matches for gear names or integers between 1 and ",
          as.character(blicc_ld$NG)
        )
      )
    } else {
      if (is.character(Gear[1])) {
        if (Gear[1] == "All")
          Gear <- 1:blicc_ld$NG
        else {
          Gear <- match(Gear, blicc_ld$gname)
          if (any(is.na(Gear))) {
            stop(
              paste0(
                "Error: Gears must be specified as 'All', or exact matches for gear names or integers between 1 and ",
                as.character(blicc_ld$NG)
              )
            )
          }
        }
      }
      Gear <- as.integer(unique(Gear))
      if (!all(dplyr::between(Gear, 1, blicc_ld$NG))) {
        stop(paste0(
          "Error: Specified gears must be between 1 and ",
          as.character(blicc_ld$NG)
        ))
      }
      return(Gear)
    }
  }
}

#' Parse the time period parameter, converting to an integer
#'
#' The time_period parameter is converted to an integer index of the periods. Exact period
#' names are also accepted.
#'
#' @inheritParams blicc_mpd
#' @param Period  A vector of 1 or more numbers or strings representing valid
#'   periods in the model
#' @return Integer period index (or stops with error message)
#' @noRd
#' 
parse_period <- function(Period, blicc_ld) {
  if (blicc_ld$NT == 1) {
    return(1L)
  } else {
    if (any(is.na(Period))) {
      stop(
        paste0(
          "Error: Periods must be specified as 'All', or exact matches for period names or integers between 1 and ",
          as.character(blicc_ld$NT)
        )
      )
    } else {
      if (is.character(Period[1])) {
        if (Period[1] == "All")
          Period <- 1:blicc_ld$NT
        else {
          Period <- match(Period, blicc_ld$tpname)
          if (any(is.na(Period))) {
            stop(
              paste0(
                "Error: Periods must be specified as 'All', or exact matches for period names or integers between 1 and ",
                as.character(blicc_ld$NT)
              )
            )
          }
        }
      }
      Period <- as.integer(unique(Period))
      if (!all(dplyr::between(Period, 1, blicc_ld$NT))) {
        stop(paste0(
          "Error: Specified periods must be between 1 and ",
          as.character(blicc_ld$NT)
        ))
      }
      return(Period)
    }
  }
}


#' Get indices for all referenced selectivity functions for the specified gears
#'
#' The gear parameter is assumed to be an integer vector already parsed. The
#' function simply returns unique indices for every selectivity function that is
#' referenced by `gear`. This allows for mixtures.
#'
#' @inheritParams parse_gear
#' @return Integer selectivity function index
#' @noRd
#'
get_selectivities <- function(gear, blicc_ld)  {
  sindx <- blicc_ld$GSbase[gear]
  mi <- (gear-1L)*2L + 1L
  mix <- mi[blicc_ld$GSmix1[mi]>0]
  for (mii in mix) {
    sindx <- with(blicc_ld, c(sindx, GSmix2[GSmix1[mii]:GSmix1[mii+1L]]))
  }
  return(sort(unique(sindx))) # indices of selectivity functions being used by Gear
}



#' Get selectivity parameter names
#'
#' The selectivity parameters listed in the data object are given names based on
#' the selectivity function number of the purpose of the parameter (location or
#' slope).
#'
#' @inheritParams parse_gear
#' @return Vector of selectivity parameter names
#' @noRd
#' 
get_sel_par_names <- function(blicc_ld)  {
  par_names <- character(blicc_ld$NP)
  j <- 1
  for (i in 1:blicc_ld$NS) {
    if (blicc_ld$fSel[i] %in% 1:3) { 
      par_names[j] <- paste0("S", as.character(i), "_loc")
      j <- j + 1
      par_names[j] <- paste0("S", as.character(i), "_slp")
      j <- j + 1
    } else {
      par_names[j] <- paste0("S", as.character(i), "_loc")
      j <- j + 1
      par_names[j] <- paste0("S", as.character(i), "_slp1")
      j <- j + 1
      par_names[j] <- paste0("S", as.character(i), "_slp2")
      j <- j + 1
    }
  }
  return(par_names)
}



