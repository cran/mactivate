\name{f_mactivate}
\alias{f_mactivate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Map Activation Layer and Inputs to Polynomial Model Inputs
}
\description{
Passes activation inputs, \code{U} into activation layer, \code{W}, to obtain new polynomial model inputs.
}
\usage{
f_mactivate(U, W)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{U}{
Numeric matrix, \code{N} x \code{d_u} of activation inputs.
}
  \item{W}{
Numeric matrix, \code{d_u} x \code{m}, the multiplicative activation layer.
}
}
\details{
This function calculates the multiplicative activations; it maps selected inputs, \code{U}, back into the input space using the m-activation layer(s).  In practice, the arg \code{W}, will be a fitted value, as created by the fitting functions.
}
\value{
Numeric matrix, \code{N} x \code{m}.  Referred to as \code{Xstar} elsewhere in this documentation.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

\donttest{

library(mactivate)

set.seed(777)


d <- 7
N <- 15000

X <- matrix(rnorm(N*d, 0, 1), N, d) ####

colnames(X) <- paste0("x", I(1:d))

############# primary effects
b <- rep_len( c(-1/4, 1/4), d )



###########

xxA <- (X[ , 1]+1/3) * (X[ , 1]-1/3) * (X[ , 3]+1/3)
xxB <- (X[ , 2]+0) * (X[ , 2]+1/3) * (X[ , 3]-0) * (X[ , 3]-1/3)
xxC <- (X[ , 3]+1/3) * (X[ , 3]-1/3)

ystar <-
X \%*\% b +
1/3 * xxA -
1/2 * xxB +
1/3 * xxC


#############

xs2 <- "y ~ . "

xtrue_formula <- eval(parse(text=xs2))

xnoint_formula <- eval(parse(text="y ~ . - xxA - xxB - xxC"))



yerrs <- rnorm(N, 0, 3)

y <- ystar + yerrs

########## standardize X
Xall <- t( ( t(X) - apply(X, 2, mean) ) / apply(X, 2, sd) )
yall <- y
Nall <- N


####### fold index
xxfoldNumber <- rep_len(1:2, N)

ufolds <- sort(unique(xxfoldNumber)) ; ufolds


############### predict
############### predict


dfx <- data.frame("y"=yall, Xall, xxA, xxB, xxC)

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
###### takes about 2 minutes

xcmact_hybrid <-
f_control_mactivate(
param_sensitivity = 10^12,
bool_free_w       = TRUE,
w0_seed           = 0.1,
w_col_search      = "alternate",
max_internal_iter = 500, #####
ss_stop           = 10^(-14), ###
escape_rate       = 1.005,
Wadj              = 1/1,
force_tries       = 0,
lambda            = 0/10000, ###
tol               = 10^(-14) ###
)




#### Fit

m_tot <- 7

Uall <- cbind(Xall, Xall)
colnames(Uall) <- paste0(rep(c("a_", "b_"), each=d), colnames(Uall))

head(Uall)

xthis_fold <- ufolds[ 1 ]


xndx_test <- which( xxfoldNumber \%in\% xthis_fold )
xndx_train <- setdiff( 1:Nall, xndx_test )

X_train <- Xall[ xndx_train, , drop=FALSE ]
y_train <- yall[ xndx_train ]
U_train <- Uall[ xndx_train, , drop=FALSE ]

xxnow <- Sys.time()
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
cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )



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

xthis_fold <- ufolds[ 1 ]

xndx_test <- which( xxfoldNumber \%in\% xthis_fold )
xndx_train <- setdiff( 1:Nall, xndx_test )

xlm <- lm(xtrue_formula , data=dfx[ xndx_train, ])
yhat <- predict(xlm, newdata=dfx[ xndx_test, ])

sqrt( mean( (y_test - yhat)^2 ) )


################ hatXstar

X_test <- Xall[ xndx_test, ]
y_test <- yall[ xndx_test ]

Xstar_test <- f_mactivate(U=U_test, W=xxls_out[[ length(xxls_out) ]][[ "What" ]])
Xi <- cbind(X_test, Xstar_test)
xlm <- lm(y_test ~ Xi)

sumxlm <- summary(xlm)
print(sumxlm)

xcoefs <- sumxlm$coefficients
xcoefs <- xcoefs[ (2+d):nrow(xcoefs), ] ; xcoefs

xndox_cu <- which( abs(xcoefs[ , "t value"]) > 3 ) ; xndox_cu


bWhat <- xxls_out[[ length(xxls_out) ]][[ "What" ]][ ,  xndox_cu ]
bWhat

wwmag <- apply(bWhat, 1, function(x) { return(sum(abs(x)))} ) ; wwmag

plot(wwmag, type="h", lwd=4,
ylim=c(0, max(wwmag)),
main="W Coefficient Total Magnitute vs Input Term",
xlab="Column of U",
ylab="Sum of magnitudes in fitted W",
cex.lab=1.3
)



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
