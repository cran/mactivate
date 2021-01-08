\name{predict.mactivate_fit_hybrid_01}
\alias{predict.mactivate_fit_hybrid_01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict from Fitted Hybrid Model
}
\description{
Predict using fitted model returned by \code{\link{f_fit_hybrid_01}}.
}
\usage{
\method{predict}{mactivate_fit_hybrid_01}(object, X0, U0=NULL, mcols, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A list of class 'mactivate_fit_hybrid_01' as returned by \code{\link{f_fit_hybrid_01}}.
}
  \item{X0}{
Numeric matrix, \code{N} x \code{d}.  Model `primary effect' inputs.
}
  \item{U0}{
Numeric matrix with \code{N} rows.  Inputs to pass to activation layer.
}
  \item{mcols}{
Scalar non-negative integer specifying which first columns of \code{W} to use.
}
  \item{\dots}{
Nothing else is required for this extension of the predict() function.
}
}
\details{
If \code{U0} is not provided, \code{X0} will be passed to activation layer.
}
\value{
\code{yhat}.  Numeric vector of length \code{N}.
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
\code{\link{f_fit_hybrid_01}}.
}
\examples{
####### Please see examples in the fitting functions
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
