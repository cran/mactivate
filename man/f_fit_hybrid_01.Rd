\name{f_fit_hybrid_01}
\alias{f_fit_hybrid_01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Multivariate Regression Model with mactivate Using Hybrid Method
}
\description{
Use hybrid algorithm (essentially a flavor of EM) to locate model parameters, i.e., primary effects, multiplicative effects, and activation parameters, \code{W}.
}
\usage{
f_fit_hybrid_01(
X, 
y, 
m_tot, 
U = NULL, 
m_start = 1, 
mact_control = f_control_mactivate(), 
verbosity = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Numerical matrix, \code{N} x \code{d} of model inputs.  Do not include intercept term.
}
  \item{y}{
Numerical vector of length \code{N}.  Model response, or output.
}
  \item{m_tot}{
Scalar non-negative integer.  Total number of columns of activation layer, \code{W}, over which to fit.
}
  \item{U}{
Numerical matrix, \code{N} x \code{d_u} of model inputs to send to the activation layer, \code{W}.  The default, \code{NULL}, instructs this function to simply use arg \code{X}.
}
  \item{m_start}{
Currently not used.
}
  \item{mact_control}{
Named list of class \code{control_mactivate_obj} as created by fun \code{\link{f_control_mactivate}} --- fitting hyperparameters.
}
  \item{verbosity}{
Scalar integer.
}
}
\details{
Please make sure to read Details in \code{\link{f_dmss_dW}} help page before using this function.
}
\value{
%A named list of class \code{mactivate_fit_hybrid_01}.
An unnamed list of class \code{mactivate_fit_hybrid_01} of length \code{m_tot + 1}.  Each node is a named list containing fitted parameter estimates.  The first top-level node of this object contains parameter estimates when fitting `primary effects' only (\code{W} has no columns), the second, parameter estimates for fitting with 1 column of \code{W}, and so on.
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
Essentially equivalent to, but likely faster than: \code{\link{f_fit_gradient_01}}.  See \code{\link{f_fit_gradient_logistic_01}} for logistic data (binomial response).
}
\examples{

xxnow <- Sys.time()

library(mactivate)

set.seed(777)


d <- 4
N <- 2000

X <- matrix(rnorm(N*d, 0, 1), N, d) ####

colnames(X) <- paste0("x", I(1:d))

############# primary effects
b <- rep_len( c(-1/2, 1/2), d )



###########

xxA <- (X[ , 1]+1/3) * (X[ , 2]-1/3)
#xxA <- (X[ , 1]+0/3) * (X[ , 2]-0/3)


ystar <-
X \%*\% b +
2 * xxA



#############





xs2 <- "y ~ . "


xtrue_formula <- eval(parse(text=xs2))

xnoint_formula <- eval(parse(text="y ~ . - xxA"))



yerrs <- rnorm(N, 0, 3)

y <- ystar + yerrs

## y <- (y - mean(y)) / sd(y)


########## standardize X
Xall <- t( ( t(X) - apply(X, 2, mean) ) / apply(X, 2, sd) )
yall <- y
Nall <- N


####### fold index
xxfoldNumber <- rep_len(1:2, N)

ufolds <- sort(unique(xxfoldNumber)) ; ufolds


############### predict
############### predict


dfx <- data.frame("y"=yall, Xall, xxA)

tail(dfx)



################### incorrectly fit LM: no interactions
xlm <- lm(xnoint_formula , data=dfx)
summary(xlm)
yhat <- predict(xlm, newdata=dfx)
sqrt( mean( (yall - yhat)^2 ) )



################### correctly fit LM
xlm <- lm(xtrue_formula, data=dfx)
summary(xlm)
yhat <- predict(xlm, newdata=dfx)
sqrt( mean( (yall - yhat)^2 ) )





################ fit using hybrid m-activation
######

m_tot <- 4

xcmact_hybrid <-
f_control_mactivate(
param_sensitivity = 10^12,
bool_free_w       = TRUE,
w0_seed           = 0.1, ### 0.01
w_col_search      = "alternate",
max_internal_iter = 500, #####
ss_stop           = 10^(-11), ###
escape_rate       = 1.02, ### 1.05
Wadj              = 1/1,
force_tries       = 0,
lambda            = 0/10000, ### hybrid only
tol               = 10^(-11) ### hybrid only
)

#### Fit

Uall <- Xall

head(Uall)

xthis_fold <- ufolds[ 1 ]

xndx_test <- which( xxfoldNumber \%in\% xthis_fold )
xndx_train <- setdiff( 1:Nall, xndx_test )

X_train <- Xall[ xndx_train, , drop=FALSE ]
y_train <- yall[ xndx_train ]
U_train <- Uall[ xndx_train, , drop=FALSE ]


xxls_out <-
f_fit_hybrid_01(
X = X_train,
y = y_train,
m_tot = m_tot,
U = U_train,
m_start = 1,
mact_control = xcmact_hybrid,
verbosity = 1
)



######### check test error

U_test <- Uall[ xndx_test, , drop=FALSE ]
X_test <- Xall[ xndx_test, , drop=FALSE ]
y_test <- yall[ xndx_test ]


yhatTT <- matrix(NA, length(xndx_test), m_tot+1)

for(iimm in 0:m_tot) {
    yhat_fold <- predict(object=xxls_out, X0=X_test, U0=U_test, mcols=iimm )
    yhatTT[ , iimm + 1 ] <- yhat_fold
}

errs_by_m <- NULL
for(iimm in 1:ncol(yhatTT)) {
    yhatX <- yhatTT[ , iimm]
    errs_by_m[ iimm ] <- sqrt(mean( (y_test - yhatX)^2 ))
    cat(iimm, "::", errs_by_m[ iimm ])
}



plot(0:(length(errs_by_m)-1), errs_by_m, type="l", xlab="m", ylab="RMSE Cost")






##################
xtrue_formula_use <- xtrue_formula


xlm <- lm(xnoint_formula , data=dfx[ xndx_train, ])
yhat <- predict(xlm, newdata=dfx[ xndx_test, ])
cat("\n\n", "No interaction model RMSE:", sqrt( mean( (y_test - yhat)^2 ) ), "\n")


xlm <- lm(xtrue_formula_use , data=dfx[ xndx_train, ])
yhat <- predict(xlm, newdata=dfx[ xndx_test, ])
cat("\n\n", "'true' model RMSE:", sqrt( mean( (y_test - yhat)^2 ) ), "\n")


cat( "Runtime:", difftime(Sys.time(), xxnow, units="secs"), "\n" )

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
