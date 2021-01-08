\name{f_control_mactivate}
\alias{f_control_mactivate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set Fitting Hyperparameters
}
\description{
Allows user a single function to tune the mactivate fitting algorithms, \code{\link{f_fit_gradient_01}}, \code{\link{f_fit_hybrid_01}}, \code{\link{f_fit_gradient_logistic_01}}.
}
\usage{
f_control_mactivate(
param_sensitivity = 10^9, 
bool_free_w = FALSE, 
w0_seed = 0.1, 
max_internal_iter = 500, 
w_col_search = "one", 
bool_headStart = FALSE, 
ss_stop = 10^(-8), 
escape_rate = 1.004, 
step_size = 1/100, 
Wadj = 1/1, 
force_tries = 0, 
lambda = 0, 
tol = 10^(-8))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param_sensitivity}{
Large positive scalar numeric.
}
  \item{bool_free_w}{
Scalar logical.  Allow values of \code{W} to wander outside [0,1]?
}
  \item{w0_seed}{
Scalar numeric.  Usually in [0,1].  Initial value(s) for multiplicative activation layer, \code{W}.
}
  \item{max_internal_iter}{
Scalar non-negative integer. \bold{Hybrid only}.  How many activation descent passes to make before refitting primary effects.
}
  \item{w_col_search}{
Scalar character.  When \code{one}, locating \code{W} and corresponding coefficients is done (progressively) one column at a time; when \code{all}, locating \code{W} and corresponding coefficients is done for current column and all previous columns; When \code{alternate}, locating \code{W} and corresponding coefficients is done (progressively) one column at a time, however, after each column is fitted, an additonal pass is made fitting current column and all previous columns.
}
  \item{bool_headStart}{
Scalar logical. \bold{Gradient only}. When \code{TRUE}, fitting first locates initial primary effects as a ``head start'' to the subsequent gradient fitting.
}
  \item{ss_stop}{
Small positive scalar numeric.  Convergence tolerance.
}
  \item{escape_rate}{
Scalar numeric no less than one and likely no greater than, say, 1.01.  Affinity for exiting a column search over \code{W}.  E.g., if 1, fitting may take a long time.  If 1.01, search for each column \code{W} will terminate relatively quickly.
}
  \item{step_size}{
Positive scalar numeric.  Initial gradient step size (in both gradient and hybrid fitting algorithms) for all parameters.
}
  \item{Wadj}{
Positive scalar numeric.  Control gradient step size (in both gradient and hybrid fitting algorithms) of \code{W}.
}
  \item{force_tries}{
Scalar non-negative integer.  Force a minimum number of fitting recursions.
}
  \item{lambda}{
Scalar numeric.  Ridge regularizer.  The actual diagonal loading imposed upon the precision matrix is equal to \code{lambda} times its original diagonal.  A value of \code{0} applies no loading; a value of \code{1} doubles the diagonal values of the precision matrix.  This is applied to primary effects only.  With gradient MLR fitting, i.e., \code{\link{f_fit_gradient_01}}, this only applies when arg \code{bool_headStart} is set to \code{TRUE} (otherwise there'd be nothing to regularize).  With hybrid MLR fitting, i.e., \code{\link{f_fit_hybrid_01}}, this regularization is applied at each LS step (see About vignette).  With logistic fitting, this arg \bold{does nothing}.  Note that with logistic fitting, we can always add a small amount of white noise to \code{X}.
}
  \item{tol}{
Small positive scalar numeric. \bold{Hybrid only}. Similar to arg \code{ss_stop} above, but controls convergence tolerance after both recursions in hybrid fitting have completed.
}
}
\details{
Fitting a mactivate model to data can/will be dramatically affected by these tuning hyperparameters.  
On one extreme, one set of hyperparameters may result in the fitting algorithm fruitlessly exiting almost immediately.  Another set of hyperparameters may send the fitting algorithm to run and run for hours.
While an ideal hyperparameterization will expeditiously fit the data.
}
\value{
Named list to be passed to \code{mact_control} arg in fitting functions.
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
\code{\link{f_fit_gradient_01}}, \code{\link{f_fit_hybrid_01}}, \code{\link{f_fit_gradient_logistic_01}}.
}
\examples{

\donttest{

library(mactivate)

set.seed(777)

d <- 20
N <- 50000

X <- matrix(rnorm(N*d, 0, 1), N, d)

colnames(X) <- paste0("x", I(1:d))

############# primary effect slopes
b <- rep_len( c(-1, 1), d )


ystar <-
X \%*\% b +
1 * (X[ , 1]) * (X[ , 2]) * (X[ , 3]) -
1 * (X[ , 2]) * (X[ , 3]) * (X[ , 4]) * (X[ , 5])

Xall <- X

errs <- rnorm(N, 0, 1)
errs <- 3 * (errs - mean(errs)) / sd(errs)

sd(errs)

y <- ystar + errs ### response

yall <- y
Nall <- N



############# hybrid example


### this control setting will exit too quickly
### compare this with example below

xcmact <-
f_control_mactivate(
param_sensitivity = 10^5,
w0_seed           = 0.1,
max_internal_iter = 500,
w_col_search      = "one",
ss_stop           = 10^(-5),
escape_rate       = 1.01,
Wadj              = 1/1,
lambda            = 1/1000,
tol               = 10^(-5)
)


m_tot <- 4

Uall <- Xall

xxnow <- Sys.time()

xxls_out <-
f_fit_hybrid_01(
X = Xall,
y = yall,
m_tot = m_tot,
U = Uall,
m_start = 1,
mact_control = xcmact,
verbosity = 1
)

cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )

yhatG <- predict(object=xxls_out, X0=Xall, U0=Uall, mcols=m_tot )

sqrt( mean( (yall  -  yhatG)^2 ) )





####################### this control setting should fit
####################### (will take a few minutes)

xcmact <-
f_control_mactivate(
param_sensitivity = 10^10, ### make more sensitive
w0_seed           = 0.1,
max_internal_iter = 500,
w_col_search      = "one",
ss_stop           = 10^(-14), ### make stopping insensitive
escape_rate       = 1.001, #### discourage quitting descent
Wadj              = 1/1,
lambda            = 1/10000,
tol               = 10^(-14) ### make tolerance very small
)


m_tot <- 4

Uall <- Xall

xxnow <- Sys.time()

xxls_out <-
f_fit_hybrid_01(
X = Xall,
y = yall,
m_tot = m_tot,
U = Uall,
m_start = 1,
mact_control = xcmact,
verbosity = 1
)

cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )

yhatG <- predict(object=xxls_out, X0=Xall, U0=Uall, mcols=m_tot )

sqrt( mean( (yall  -  yhatG)^2 ) )


xxls_out

Xstar <- f_mactivate(U=Uall, W=xxls_out[[ m_tot+1 ]][[ "What" ]])
colnames(Xstar) <- paste0("xstar_", seq(1, m_tot))
Xall <- cbind(Xall, Xstar)

xlm <- lm(yall~Xall)
summary(xlm)

}

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
