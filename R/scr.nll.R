#=====================================#
#       Log-likelihood function       #
#=====================================#
#' @export
scr.nll = function(pars, caps, traps, mask, maskDists) {
  scr_nll(pars, caps, traps, mask, maskDists)
}

#==========================================================================#
#==========================================================================#

#=====================================#
#  Log-likelihood function (acoustic) #
#=====================================#
#' @export
scr.nll.acoustic = function(pars, caps, traps, mask, maskDists) {
  scr_nll_acoustic(pars,
                   caps[, -ncol(caps)],
                   traps,
                   mask,
                   maskDists,
                   table(caps[, ncol(caps)]))
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
