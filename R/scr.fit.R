#=====================================#
#       Main fitting function         #
#=====================================#
#' @export
scr.fit = function(capthist, traps, mask,
                   start = NULL, acoustic = FALSE) {
  ## General error/exception handling
  if(length(start) == 4 && acoustic == FALSE) {
    warning("Data treated as acoustic captures")
  }

  ## Checking to see if things need unpacking
  if(class(capthist) == "capthist") {
    capthist = capthist[, 1, ]
    traps = capthist$traps
  }

  ## Transforming the start values
  ## - Note that acoustic captures have one extra parameter (lambda_c)
  if(acoustic) {
    start = c(log(start[1]),
            qlogis(start[2]),
            log(start[3]),
            start[4])
  } else {
    start = c(log(start[1]),
              qlogis(start[2]),
              log(start[3]))
  }

  ## Calculating mask distances before giving to optim
  ## - More efficient
  maskDists = eucdist_nll(mask, traps)

  ## Taking the capthist, traps, and mask and maximising likelihood
  ## - Note that likelihood function changes for acoustic
  ## - Requires start values
  if(acoustic) {
    fit = optim(start, scr.nll.acoustic,
                caps = capthist,
                traps = traps,
                mask = mask,
                maskDists = maskDists)
  } else {
    fit = optim(start, scr.nll,
                caps = capthist,
                traps = traps,
                mask = mask,
                maskDists = maskDists)
  }

  ## Returning the fitted parameters in a named vector
  fittedPars = c(fit$par)
  parNames = ifelse(acoustic,
                    c("D", "g0", "sigma"),
                    c("D", "lambda_0", "sigma"))
  setNames(c(exp(fittedPars[1]),
             plogis(fittedPars[2]),
             exp(fittedPars[3])),
           nm = parNames)
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
