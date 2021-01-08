### R code from vignette source 'mactivation_tutorial_01.Snw'

###################################################
### code chunk number 1: mactivation_tutorial_01.Snw:108-111
###################################################
options(width=49)
options(prompt=" ")
options(continue="   ")


###################################################
### code chunk number 2: a1
###################################################

library(mactivate)
set.seed(777)

## tiny
d <- 11
N <- 3000

X <- matrix(rnorm(N*d, 1, 1), N, d)
colnames(X) <- paste0("x", I(1:d))

b <- rep_len( c(-1, 1), d )

ystar <-
X %*% b +
1/3 * X[ , 1] * X[ , 2] * X[ , 3] -
1/3 * X[ , 3] * X[ , 4] * X[ , 5] * X[ , 6] +
1/2 * X[ , 8] * X[ , 9] -
2   * X[ , 1] * X[ , 2] * X[ , 7] * X[ , 11]

xtrue_formula <- eval(parse(text="y ~ . + x1:x2:x3 + x3:x4:x5:x6 + x8:x9 + x1:x2:x7:x11"))

xnoint_formula <- eval(parse(text="y ~ ."))

errs <- rnorm(N, 0, 3)

y <- ystar + errs

Xall <- X
yall <- y
Nall <- N

dfx <- data.frame("y"=yall, Xall)


###################################################
### code chunk number 3: a2
###################################################
xlm <- lm(y ~ . , data=dfx)
yhat <- predict(xlm, newdata=dfx)
sqrt( mean( (yall -  yhat)^2 ) )


###################################################
### code chunk number 4: a3
###################################################
xlm <- lm(xtrue_formula , data=dfx)
yhat <- predict(xlm, newdata=dfx)
sqrt( mean( (yall  -  yhat)^2 ) )


###################################################
### code chunk number 5: a4
###################################################
xcmact_hybrid <-
f_control_mactivate(
param_sensitivity = 10^10,
w0_seed           = 0.1,
w_col_search      = "one",
bool_headStart    = FALSE, ### gradient
max_internal_iter = 500, ##### small -- exits automatically, don't set this too small
ss_stop           = 10^(-8), ### small
escape_rate       = 1.01,
Wadj              = 1/1,
tol               = 10^(-8)
)


###################################################
### code chunk number 6: a5
###################################################

m_tot <- 5

Uall <- Xall

xxnow <- Sys.time()

xxls_out <-
f_fit_hybrid_01(
X = Xall,
y = yall,
m_tot = m_tot,
U = Uall,
m_start = 1,
mact_control = xcmact_hybrid,
verbosity = 5
)

cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )


###################################################
### code chunk number 7: a6
###################################################
class(xxls_out)
yhatall <- predict(object=xxls_out, X0=Xall, U0=Uall, mcols=m_tot)
sqrt( mean( (yall - yhatall)^2 ) )


