\name{scr.fit}
\alias{scr.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial Capture-Recapture
}
\description{
Maximises likelihood and returns a set of parameter estimates, given a matrix of capture histories, trap locations, and mask points.
}
\usage{
scr.fit(capthist, traps, mask, start = NULL, acoustic = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{capthist}{
	Matrix of capture histories.
  }
  \item{traps}{
  Matrix of trap location coordinates.
  }
  \item{mask}{
  Matrix of mask point coordinates.
  }
  \item{start}{
  Vector of starting values for the optimisation.
  }
  \item{acoustic}{
  Fit acoustic capture data. See 'Details'.
  }
  \item{binom}{
  Specifies whether data are binomial. See 'Details'.
  }
  \item{toa}{
  Time of arrival matrix. See 'Details'.
  }
  \item{speed_sound}{
  Speed of sound. Default is 330(m/s).
  }
}
\details{
Parameter estimates are obtained by maximising likelihood (via  \code{optim}).

When fitting cue-based (e.g. acoustic) data, \code{capthist} matrices must contain column of 'animal IDs' as their last column.

Fitting binary data means that the parameter \code{g0} is included, rather than \code{lambda0}. Starting parameter values will thus default to different values, depending on whether the data are binary.

Time of arrival matrix has the same dimensions as the capture history matrix. Each entry is now the time elapsed between the beginning of recording, and a detection. Inclusion of this matrix requires an extra parameter appended to \code{start}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Borchers, D. L. and Efford, M. G. (2008) Spatially explicit maximum likelihood methods for capture--recapture studies. \emph{Biometrics} \strong{64}, 377--385.
}
\author{
Callum Young
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[secr]{secr}},
\code{\link[scr]{scr.sim}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.

## Loading in the example data (including trap and mask locations).
data(scr.data)

## Then running scr.sim() with lambda_0 = 5, sigma = 75, density = 2.
capt = scr.sim(5, 75, scr.data$traps, density = 2)

## Starting parameters are saved for use in the model fitting function.
par.start = c(1, 1, 50)

## Finally, fitting a model to the data,
##  using the trap and mask points in the example data.
scr.fit(capt, scr.data$traps, scr.data$mask, start = par.start)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
