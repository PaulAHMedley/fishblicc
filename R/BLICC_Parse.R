# BLICC Parsing Functions ---------------------------------------------------

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
#' The selectivity function parameter is converted to an integer index of the
#' function, or `NA` is returned with an error message.
#'
#' @inheritParams blicc_mpd
#' @param sel_fun  A number or string representing a valid available function.
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

    if (sel_fun[1]=="all")
      return(1L:blicc_ld$NG)

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



#' Get indices for all referenced selectivity functions for the specified gears
#'
#' The gear parameter is assumed to be an integer vector already parsed. The
#' function simply returns unique indices for every selectivity function that is
#' referenced by Gear. This allows for mixtures.
#'
#' @inheritParams parse_gear
#' @return Integer selectivity function index
#' @noRd
#'
get_selectivities <- function(Gear, blicc_ld)  {
  sindx <- blicc_ld$GSbase[Gear]
  mi <- (Gear-1L)*2L + 1L
  mix <- mi[blicc_ld$GSmix1[mi]>0]
  for (mii in mix) {
    sindx <- with(blicc_ld, c(sindx, GSmix2[GSmix1[mii]:GSmix1[mii+1L]]))
  }
  return(sort(unique(sindx))) # indices of selectivity functions being used by Gear
}


