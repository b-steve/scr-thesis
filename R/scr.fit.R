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
  ## - Machine minimum subtracted so as to avoid log errors
  start = start - .Machine$double.xmin
  if(acoustic) {
    start = c(log(start[1]),
            qlogis(start[2]),
            log(start[3]),
            log(start[4]))
  } else {
    ## 3 parameters, all logged
    start = log(start)
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
                maskDists = maskDists,
                hessian = TRUE)
  } else {
    fit = optim(start, scr.nll,
                caps = capthist,
                traps = traps,
                mask = mask,
                maskDists = maskDists,
                hessian = TRUE)
  }

  ## Returning the fitted parameters in a named vector
  fittedPars = fit$par
  if(acoustic) {
    parNames = c("D", "g0", "sigma", "lambda_c")

    fittedPars = c(exp(fittedPars[1]),
                   plogis(fittedPars[2]),
                   exp(fittedPars[3]),
                   exp(fittedPars[4]))
  } else {
    parNames = c("D", "lambda_0", "sigma")

    fittedPars = c(exp(fittedPars[1]),
                   exp(fittedPars[2]),
                   exp(fittedPars[3]))
  }

  ## Calculating confidence intervals
  ## - Using the (sqrt of) diagonals of (-ve) Hessian obtained from optim
  ##    - i.e. Information matrix
  ## - Wald CIs calculated by sapply() loop
  ##    - Loops through each of fitted parameters and calculates lower/upper bounds
  se = sqrt(diag(solve(fit$hess)))
  waldCI = t(sapply(1:length(fittedPars),
                    function(i) fittedPars[i] + (c(-1, 1) * (qnorm(0.975) * se[i]))))
  waldCI = t(cbind(exp(waldCI[1, ]),
                   plogis(waldCI[2, ]),
                   exp(waldCI[3:nrow(waldCI), ])))

  results = cbind(fittedPars, waldCI)
  dimnames(results) = list(parNames, c("Estimate", "Lower", "Upper"))
  results
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
