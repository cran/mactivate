\name{f_fit_gradient_logistic_01}
\alias{f_fit_gradient_logistic_01}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit Logistic Multivariate Regression Model with mactivate Using Gradient Descent
}
\description{
Use simple gradient descent to locate logistic model parameters, i.e., primary effects, multiplicative effects, and activation parameters, \code{W}.
}
\usage{
f_fit_gradient_logistic_01(
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
Integer vector of length \code{N}, elements in \code{{0, 1}}.  Binomial model response, or output.
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
%A named list of class \code{mactivate_fit_gradient_logistic_01}.
An unnamed list of class \code{mactivate_fit_gradient_logistic_01} of length \code{m_tot + 1}.  Each node is a named list containing fitted parameter estimates.  The first top-level node of this object contains parameter estimates when fitting `primary effects' only (\code{W} has no columns), the second, parameter estimates for fitting with 1 column of \code{W}, and so on.
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
See \code{\link{f_fit_gradient_01}} or \code{\link{f_fit_gradient_logistic_01}} for MLR data (numerical response).
}
\examples{

xxnow <- Sys.time()

library(mactivate)

set.seed(777)


d <- 4
N <- 2400

X <- matrix(rnorm(N*d, 0, 1), N, d) ####

colnames(X) <- paste0("x", I(1:d))

############# primary effects
b <- rep_len( c(-1/2, 1/2), d )


xxA <- (X[ , 1]+1/3) * (X[ , 2]-1/3)
xxB <- (X[ , 1]+0/3) * (X[ , 2]-0/3) * (X[ , 4]-1/3)


ystar <-
X \%*\% b +
2 * xxA -
1 * xxB

xs2 <- "y ~ . "

xtrue_formula <- eval(parse(text=xs2))

xnoint_formula <- eval(parse(text="y ~ . - xxA - xxB"))


ysigmoid <- 1 / (1 + exp(-ystar))

range(ysigmoid)

y <- rbinom(size=1, n=N ,prob=ysigmoid)


Nall <- N

cov(X)

yall <- y

Xall <- X

### Xall <- X + rnorm(prod(dim(X)), 0, 1/10000) ### add a little noise -- optional

sd(y)


dfx <- data.frame("y"=yall, Xall, xxA, xxB)

tail(dfx)



################### incorrectly fit LM: no interactions
xglm <- glm(xnoint_formula , data=dfx, family=binomial(link="logit"))
summary(xglm)
yhat <- predict(xglm, newdata=dfx, type="response")
mean( f_logit_cost(y=yall, yhat=yhat) )


####### known true
xglm <- glm(xtrue_formula , data=dfx, family=binomial(link="logit"))
summary(xglm)
yhat <- predict(xglm, newdata=dfx, type="response")
mean( f_logit_cost(y=yall, yhat=yhat) )



xxfoldNumber <- rep_len( 1:4, Nall )

ufolds <- sort(unique(xxfoldNumber))

######################

xthis_fold <- ufolds[ 1 ]

xndx_test <- which( xxfoldNumber \%in\% xthis_fold )

xndx_train <- setdiff( 1:Nall, xndx_test )

##################

X_train <- Xall[ xndx_train, , drop=FALSE ]

X_test <- Xall[ xndx_test, , drop=FALSE ]

y_train <- yall[ xndx_train ]

y_test <- yall[ xndx_test ]



###################

m_tot <- 4

xcmact_gradient <-
f_control_mactivate(
param_sensitivity = 10^11,
bool_free_w       = FALSE,
w0_seed           = 0.05,
#w_col_search      = "alternate",
w_col_search      = "one",
bool_headStart    = TRUE,
ss_stop           = 10^(-12), ### very small
escape_rate       = 1.02,
step_size         = 1,
Wadj              = 1/1,
force_tries       = 0,
lambda            = 1/1 #### does nothing here
)


Uall <- Xall


X_train <- Xall[ xndx_train, , drop=FALSE ]
y_train <- yall[ xndx_train ]


xxls_out <-
f_fit_gradient_logistic_01(
X = X_train,
y = y_train,
m_tot = m_tot,
U = X_train,
m_start = 1,
mact_control = xcmact_gradient,
verbosity = 0
)


######### check test error

U_test <- Xall[ xndx_test, , drop=FALSE ]
X_test <- Xall[ xndx_test, , drop=FALSE ]
y_test <- yall[ xndx_test ]


yhatTT <- matrix(NA, length(xndx_test), m_tot+1)

for(iimm in 0:m_tot) {
    yhat_fold <- predict(object=xxls_out, X0=X_test, U0=U_test, mcols=iimm )
    yhatTT[ , iimm + 1 ] <- yhat_fold[[ "p0hat" ]]
}

errs_by_m <- NULL
for(iimm in 1:ncol(yhatTT)) {
    yhatX <- yhatTT[ , iimm]
    errs_by_m[ iimm ] <- mean( f_logit_cost(y=y_test, yhat=yhatX) )
    cat(iimm, "::", errs_by_m[ iimm ])
}

##### plot test Logit vs m

plot(0:(length(errs_by_m)-1), errs_by_m, type="l", xlab="m", ylab="Logit Cost")





################## test off 'correct' model
xtrue_formula_use <- xtrue_formula



xglm <- glm(xnoint_formula , data=dfx[ xndx_train, ], family=binomial(link="logit"))
yhat <- predict(xglm, newdata=dfx[ xndx_test, ], type="response")
cat("\n\n", "No interaction model logit:", mean( f_logit_cost(y=y_test, yhat=yhat) ), "\n")


xglm <- glm(xtrue_formula_use , data=dfx[ xndx_train, ], family=binomial(link="logit"))
yhat <- predict(xglm, newdata=dfx[ xndx_test, ], type="response")
cat("\n\n", "'true' model logit:", mean( f_logit_cost(y=y_test, yhat=yhat) ) , "\n")


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
