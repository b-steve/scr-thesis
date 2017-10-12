#=====================================#
#        Main fitting function        #
#=====================================#
#' @export
scr.fit = function(capthist, traps, mask,
                   start = NULL) {
  ## Checking to see if things need unpacking
  if(class(capthist) == "capthist") {
    capthist = capthist[, 1, ]
    traps = capthist$traps
  }

  ## Calculating mask distances before giving to optim
  ## - More efficient
  maskDists = eucdist_nll(mask, traps)

  ## Taking the capthist, traps, and mask and maximising likelihood
  ## - Requires start values
  fit = optim(start, scr.nll,
              caps = capthist,
              traps = traps,
              mask = mask,
              maskDists = maskDists)

  ## Returning the fitted parameters in a named vector
  setNames(fit$par, nm = c("D", "g0", "sigma"))
}

#==========================================================================#
#==========================================================================#

#=====================================#
#      Acoustic fitting function      #
#=====================================#
#' @export
acoustic.scr.fit = function(capthist, traps, mask,
                            start = NULL) {
  ## Transforming start parameters
  start = c(log(start[1]),
            qlogis(start[2]),
            log(start[3]),
            start[4])

  ## Checking to see if things need unpacking
  if(class(capthist) == "capthist") {
    capthist = capthist[, 1, ]
    traps = capthist$traps
  }

  ## Calculating mask distances before giving to optim
  ## - More efficient
  maskDists = eucdist_nll(mask, traps)

  ## Taking the capthist, traps, and mask and maximising likelihood
  ## - Requires start values
  fit = optim(start, scr.nll,
              caps = capthist,
              traps = traps,
              mask = mask,
              maskDists = maskDists,
              nCalls = table(capthist[, ncol(capthist)]))

  ## Returning the fitted parameters in a named vector
  setNames(fit$par, nm = c("D", "g0", "sigma"))
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
