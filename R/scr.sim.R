#=====================================#
#      Main simulation function       #
#=====================================#
#' @export
scr.sim = function(lambda_0, sigma, traplocs,
                   density = 50,
                   distr = "pois",
                   limits = list(xlim = NULL, ylim = NULL),
                   counts = "counts",
                   draw = FALSE,
                   ...) {
  ## Setting up the total survey area
  ##  - Survey area (vs. trap area) is based on extreme trap co-ordinates
  ##    - Ranges are taken from each x and y column
  ##    - Co-ordinates are extended by (5 * sigma); stored as survey area box
  ##  - Then some error handling (i.e. only x or only y coordinates missing)
  if(is.null(limits$xlim) & is.null(limits$ylim)) {
    limits = list(xlim = range(traplocs[, 1], na.rm = TRUE),
                  ylim = range(traplocs[, 2], na.rm = TRUE))
    limits = lapply(limits, "+", c(-1, 1) * (5 * sigma))
  } else if(is.null(limits$xlim) & !is.null(limits$ylim)) {
    stop("X co-ordinates missing")
  } else if(!is.null(limits$xlim) & is.null(limits$ylim)) {
    stop("Y coordinates missing")
  }

  ## Generating random points on the area
  ##  - Area calculated from xlim and ylim
  ##  - Total number of animals N ~ Pois(DA)
  ##    - Density (D)is given per hectare, so divide by 10,000 to get per metre
  area = (limits$xlim[2] - limits$xlim[1]) * (limits$ylim[2] - limits$ylim[1])
  n = rpois(1, (density / 10000) * area)
  coords = pointgen(n, xlim = limits$xlim, ylim = limits$ylim)

  if(draw) {
    ## Setting up the survey area
    plot.new()
    plot.window(xlim = limits$xlim,
                ylim = limits$ylim,
                xaxs = "i", yaxs = "i")
    box()

    ## Setting up the traps
    points(traplocs, pch = 3, col = "red")

    ## Plotting the activity centres
    points(coords, ...)
  }

  ## Setting up the random count generation - depending on the distribution
  if(distr == "pois") {
    rDistr = paste0("r", distr, "(length(d), lambda_0 * exp(-d^2 / (2 * sigma^2)))")
  } else if(distr == "bernoulli" | distr == "binom") {
    rDistr = paste0("r", "binom", "(length(d), 1, lambda_0 * exp(-d^2 / (2 * sigma^2)))")
  } else if(distr == "negbin" | distr == "nbinom") {
    ifelse(!is.null(list(...)$size), size <- list(...)$size, size <- 2)
    rDistr = paste0("r", "nbinom", "(length(d), mu = lambda_0 * exp(-d^2 / (2 * sigma^2)), size = size)")
  }

  ## Filling the omega matrix row-by-row
  ##  - Counts are randomly generated based on the specified distribution
  ##  - A row will only be added if its sum > 0; i.e. if at least 1 of the traps had a detection.
  ##  - The simulated counts are removed at the end of each loop, just to keep things tidy.
  omega = NULL
  for(i in 1:nrow(coords)) {
    d = eucdist(coords[i,], traplocs)
    simCounts = eval(parse(text = rDistr))
    if(sum(simCounts) != 0) {
      omega = rbind(omega,
                    simCounts)
    }
    rm(simCounts)
  }

  ## Removing the row names given as a result of rbind()
  rownames(omega) = NULL

  ## Converting the count data to binary, if the count type = "binary"
  if(counts == "binary") {
    omega = ifelse(omega > 0, 1, 0)
  }

  omega
}

#==========================================================================#
#==========================================================================#

#=====================================#
#    make.capthist format function    #
#=====================================#
#' @export
toCapthist = function(captures) {
  formatted = NULL
  for(rowNum in 1:nrow(captures)) {
    ## Control flow for matrices with/without names for traps
    ##  - If captures DON'T have names for traps, will expand the trap/individual combo
    ##    - Specifically, the trap given is the trap name
    ##  - If captures have names for traps, then it'll just give the trap (column) number
    if(is.null(colnames(captures))) {
      call = "rep(which(captures[rowNum, ] != 0), captures[rowNum, which(captures[rowNum, ] != 0)])"
    } else {
      call = "rep(as.numeric(names(which(captures[rowNum, ] != 0))), captures[rowNum, which(captures[rowNum, ] != 0)])"
    }
    trapNums = eval(parse(text = call))

    ## Putting together all of the columns in the final matrix
    ## Also converting it to a data frame, with column names
    formatted = rbind(formatted, matrix(c(rep(1, length(trapNums)),
                                          rep(rowNum, length(trapNums)),
                                          rep(1, length(trapNums)),
                                          trapNums),
                                        ncol = 4))
  }
  formatted = as.data.frame(formatted)
  colnames(formatted) = c("session", "ID", "occasion", "trap")
  formatted
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
