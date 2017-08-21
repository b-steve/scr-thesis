\name{toCapthist}
\alias{toCapthist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert matrix into capthist format
}
\description{
Convert a matrix returned by \code{scr.sim} into a format expected by \code{make.capthist}.
}
\usage{
toCapthist(captures)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{captures}{
	Matrix returned by \code{scr.sim}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
## Then running scr_sim() with lambda_0 = 5, sigma = 1.25
traps = cbind(rep(seq(1, 9, length = 9), 3),
             rep(seq(9, 1, length = 9), each = 9))
traps = traps[sample(1:nrow(traps), (nrow(traps) * 1) / 5), ]

captures = scr_sim(5, 1.25, traps)

toCapthist(captures)

%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
