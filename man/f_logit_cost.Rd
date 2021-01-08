\name{f_logit_cost}
\alias{f_logit_cost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Logistic Cost
}
\description{
Calculate the logistic cost of probability predictions of a dichotomous outcome.
}
\usage{
f_logit_cost(y, yhat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
Numeric vector.  The outcome vector.  Must be in \{0, 1\}.
}
  \item{yhat}{
Numeric vector.  Prediction vector.  Should be in (0, 1) -- the open unit interval.  In an inferential setting, one should probably never make a prediction of zero or one; however, values of zero or one are allowed, provided they are ``correct''.
}
}
\details{
This function is included in this library as a convenience.  
}
\value{
A numeric vector of length equal to \code{y} and \code{yhat}.  The logistic cost associated with each corresponding prediction.
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
\code{\link{f_fit_gradient_logistic_01}}, \code{\link{predict.mactivate_fit_gradient_logistic_01}}.
}
\examples{
y <- c(0, 0, 1, 1)
yhat <- rep(1/2, length(y))

mean( f_logit_cost(y=y, yhat=yhat) )
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
