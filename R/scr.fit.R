#=====================================#
#      Main simulation function       #
#=====================================#
#' @export
scr.fit = function(capthist, traps, mask,
                   start = NULL) {

  ## Taking the capthist, traps, and mask and maximising likelihood
  ## - Requires start values
  fit = optim(start, scr.nll,
              caps = capthist,
              traps = traps,
              mask = mask)

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
