\name{mactivate-package}
\alias{mactivate-package}
\alias{mactivate}
\docType{package}
\title{
m-activation
}
\description{
\packageDescription{mactivate}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{mactivate}
\packageIndices{mactivate}
Please make sure to read Details in \code{\link{f_dmss_dW}} help page before using this library.
This package allows the user to extend the usual multivariate regression solution by adding (parallel) multiplicative ``activation layers.''  These activation layers can be very useful for identifying input interactions, and, if the user wishes, transparently test the appropriateness of input transformations.  Three functions are provided for fitting data, \code{\link{f_fit_hybrid_01}} and \code{\link{f_fit_gradient_01}} for a numeric response (usual MLR), and \code{\link{f_fit_gradient_logistic_01}} for a binary response (multivariate logistic regresssion).
The user is encouraged to consult the ``About'' vignette as well as the examples available in the respective functions' documentation for details about m-activation and practical examples of implementation.
}
\author{
\packageAuthor{mactivate}

Maintainer: \packageMaintainer{mactivate}
}
%\references{
%~~ Literature or other references for background information ~~
%}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%~~ the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## please see docs for individual functions.
}
