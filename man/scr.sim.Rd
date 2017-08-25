\name{scr.sim}
\alias{scr.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Capture Data
}
\description{
Randomly generates a spatial capture-recapture data set, and returns a matrix of encounter rates for each animal at each trap (omega).
}
\usage{
scr.sim(lambda_0, sigma, traplocs, density = 50, distr = "pois",
    limits = list(xlim = NULL, ylim = NULL),
    counts = c("count", "binary"), draw = FALSE,
    acoustic = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lambda_0}{
	Baseline encounter rate (i.e. encounter rate for d = 0)
  }
  \item{sigma}{
    Rate of decay of encounter rates (as distance between trap and animal activity centre increases).
  }
  \item{traplocs}{
    Matrix of trap location coordinates.
  }
  \item{density}{
	  Point (animal) density in animals per hectare. Default value is 50.
  }
  \item{distr}{
    Distribution of count data (for each trap) to be simulated. Default is Poisson.
  }
  \item{limits}{
	  List of length 2, giving x and y coordinates ranges. Default value is \code{NULL}; see 'Details'.
  }
  \item{counts}{
    Return matrix as counts, or as binary data.
  }
  \item{draw}{
    Plot the randomly generated points and traps. Default value is \code{FALSE}.
  }
  \item{acoustic}{
    Simulate acoustic capture histories. Default value is \code{FALSE}. See 'Details'.
  }
  \item{...}{
	  Further \link[graphics::par]{graphical parameters} may also be supplied as arguments.
  }
}
\details{
Total number of points (animal activity centres) is randomly generated from a Poisson distribution: N ~ Pois(DA), where D = \code{dens}, and A = area (automatically calculated from the trap locations; \code{xlim}, \code{ylim} can be specified manually).

Bounding box is calculated from the range given by the coordinates in \code{traplocs};points are randomly generated within this area, but not all points will be detected by the traps (depends on \code{sigma}).

Acoustic capture histories contain 1 or more rows per 'call' detected by each trap. The final column of the capture history matrix identifies each animal.
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
\code{\link[secr]{secr}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Generating some random trap locations (in a grid).
## Then running scr.sim() with lambda_0 = 5, sigma = 1.25
traps = cbind(rep(seq(1, 9, length = 9), 3),
             rep(seq(9, 1, length = 9), each = 9))
traps = traps[sample(1:nrow(traps), (nrow(traps) * 1) / 5), ]

scr.sim(5, 1.25, traps)

%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
