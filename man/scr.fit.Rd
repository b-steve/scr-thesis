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
scr.fit(capthist, traps, mask, start = NULL)
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
}
\details{
Parameter estimates are obtained by maximising likelihood (via  \code{optim}).
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
