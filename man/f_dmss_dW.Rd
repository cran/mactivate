\name{f_dmss_dW}
\alias{f_dmss_dW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Derivative of Cost Function wrt W
}
\description{
Calculate the first derivative of objective function with respect to W, given data and requisite model parameter values.
}
\usage{
f_dmss_dW(U, Xstar, W, yerrs, cc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{
Numeric matrix, \code{N} x \code{d_u} of activation inputs.
}
  \item{Xstar}{
Numeric matrix, \code{N} x \code{m}.  The ``new'' polynomial inputs created by applying the activation layer, \code{W}, to \code{U}.  Accomplished computationally with fun \code{\link{f_mactivate}}.
}
  \item{W}{
Numeric matrix, \code{d_u} x \code{m}, the multiplicative activation layer.
}
  \item{yerrs}{
Numeric vector of length \code{N}.  \code{y} minus \code{yhat}.
}
  \item{cc}{
Numeric vector of length \code{m}.  Coefficients for \code{Xstar}.
}
}
\details{
There is really no need for user to call this function directly; this function is called by the fitting functions in this library.
\bold{Important.} Computationally there are (at least) two ways to solve this derivative, one is O(Nd), the other is O(Nd^2) (d is the number of columns in \code{U}).  This function uses the first, computationally less expensive method.  It is not an approximation; the simplification occurs simply by dividing out the appropriate partial term rather than taking the full product of terms across \code{U}.  This has a very important implication of which we must be aware: \bold{zeros in \code{U} may result in division by zero!}  This function will handle the errors, but the ultimate consequence of zeros in \code{U} is that the derivative returned by this function may not be accurate.  We should eliminate zeros in \code{U}.  Standardizing \code{U} is one good solution.  If zeros are only present because of ``one-hot'' indicators (dummies), another possible solution is to substitute -1 for 0 (actually not a bad practice anyway).
}
\value{
Numeric matrix, \code{d_u} x \code{m}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{f_mactivate}}
}
\examples{
#######
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
