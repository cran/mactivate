### R code from vignette source 'mactivation_examples_02.Snw'

###################################################
### code chunk number 1: mactivation_examples_02.Snw:107-110
###################################################
options(width=63)
options(prompt=" ")
options(continue="   ")


###################################################
### code chunk number 2: a1 (eval = FALSE)
###################################################
## 
## library(mactivate)
## 
## set.seed(777)
## 
## 
## d <- 56
## N <- 100000
## 
## X <- matrix(rnorm(N*d, 0, 1), N, d) ####
## 
## colnames(X) <- paste0("x", I(1:d))
## 
## #############  primary effects
## b <- rep_len( c(-1/4, 1/4), d )
## 
## 
## 
## ###########
## 
## ystar <-
## X %*% b +
## 1/3 * (X[ , 11]+1) * (X[ , 12]-1) * (X[ , 33]+1) -
## 1/2 * (X[ , 30]+0) * (X[ , 44]+1) * (X[ , 45]-0) * (X[ , 56]-1) +
## 1/3 * (X[ , 6]+1) * (X[ , 47]-1) -
## 1/2 * (X[ , 1]-1) * (X[ , 32]+0) * (X[ , 33]+1) * (X[ , 34]-0) * 
## (X[ , 45]-0) * (X[ , 51]-1)
## 
## m_tot <- 12
## #############
## 
## 
## 
## xs1 <- 
## "y ~ . + x11:x12:x33 + x30:x44:x45:x56 + x6:x47 + x1:x32:x33:x34:x45:x51"
## 
## xs2 <- 
## "y ~ . + x11*x12*x33 + x30*x44*x45*x56 + x6*x47 + x1*x32*x33*x34*x45*x51"
## 
## 
## xnotQuiteTrue_formula <- eval(parse(text=xs1))
## xtrue_formula <- eval(parse(text=xs2))
## 
## xnoint_formula <- eval(parse(text="y ~ ."))
## 
## 
## 
## yerrs <- rnorm(N, 0, 3)
## 
## y <- ystar + yerrs
## 
## ## y <- (y - mean(y)) / sd(y)
## 
## 
## ########## standardize X
## Xall <- t( ( t(X) - apply(X, 2, mean) ) / apply(X, 2, sd) )
## yall <- y
## Nall <- N
## 
## 
## ####### fold index
## xxfoldNumber <- rep_len(1:2, N)
## 
## ufolds <- sort(unique(xxfoldNumber)) ; ufolds
## 
## 
## ############### predict
## ############### predict
## 
## dfx <- data.frame("y"=yall, Xall)
## 
## 
## 
## ################### incorrectly fit LM: no interactions
## 
## xlm <- lm(xnoint_formula , data=dfx)
## summary(xlm)
## yhat <- predict(xlm, newdata=dfx)
## sqrt( mean( (yall - yhat)^2 ) )
## 
## 
## ################### incorrectly fit LM: no lower-order interactions
## xlm <- lm(xnotQuiteTrue_formula , data=dfx)
## summary(xlm)
## yhat <- predict(xlm, newdata=dfx)
## sqrt( mean( (yall - yhat)^2 ) )
## 
## 
## ################### correctly fit LM
## xlm <- lm(xtrue_formula, data=dfx)
## summary(xlm)
## yhat <- predict(xlm, newdata=dfx)
## sqrt( mean( (yall - yhat)^2 ) )
## 
## 
## 
## 
## 
## ################ fit using hybrid m-activation
## ###### takes about 1.5-3 hours
## 
## xcmact_hybrid <-
## f_control_mactivate(
## param_sensitivity = 10^12,
## bool_free_w = TRUE,
## w0_seed = 0.01,
## w_col_search = "alternate",
## max_internal_iter = 500, #####
## ss_stop = 10^(-14), ###
## escape_rate = 1.003,  #### 1.002,
## Wadj = 1/1,
## force_tries = 0,
## lambda = 0/10000, ### hybrid only
## tol = 10^(-14) ### hybrid only
## )
## 
## 
## 
## 
## #### Fit
## 
## 
## Uall <- Xall
## 
## xthis_fold <- ufolds[ 1 ]
## 
## 
## 
## xndx_test <- which( xxfoldNumber %in% xthis_fold )
## xndx_train <- setdiff( 1:Nall, xndx_test )
## 
## X_train <- Xall[ xndx_train, , drop=FALSE ]
## y_train <- yall[ xndx_train ]
## 
## xxnow <- Sys.time()
## xxls_out <-
## f_fit_hybrid_01(
## X = X_train,
## y = y_train,
## m_tot = m_tot,
## U = X_train,
## m_start = 1,
## mact_control = xcmact_hybrid,
## verbosity = 1
## )
## cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )
## 
## 
## 
## 
## 
## ######### check test error
## 
## U_test <- Xall[ xndx_test, , drop=FALSE ]
## X_test <- Xall[ xndx_test, , drop=FALSE ]
## y_test <- yall[ xndx_test ]
## 
## 
## yhatTT <- matrix(NA, length(xndx_test), m_tot+1)
## 
## for(iimm in 0:m_tot) {
##     yhat_fold <- predict(object=xxls_out, X0=X_test, U0=U_test, mcols=iimm )
##     yhatTT[ , iimm + 1 ] <- yhat_fold
## }
## 
## errs_by_m <- NULL
## for(iimm in 1:ncol(yhatTT)) {
##     yhatX <- yhatTT[ , iimm]
##     errs_by_m[ iimm ] <- sqrt(mean( (y_test - yhatX)^2 ))
##     cat(iimm, "::", errs_by_m[ iimm ])
## }
## 
## 
## ##### plot test RMSE vs m
## 
## plot(0:(length(errs_by_m)-1), errs_by_m, type="l", xlab="m", ylab="RMSE Cost")
## 
## #### 11 best
## 
## 
## 
## 
## ################## known 'true' for non zero-centered is
## xtrue_formula_use <- xtrue_formula
## 
## xthis_fold <- ufolds[ 1 ]
## 
## xndx_test <- which( xxfoldNumber %in% xthis_fold )
## xndx_train <- setdiff( 1:Nall, xndx_test )
## 
## xlm <- lm(xtrue_formula_use , data=dfx[ xndx_train, ])
## yhat <- predict(xlm, newdata=dfx[ xndx_test, ])
## 
## 
## 
## 
## sqrt( mean( (y_test - yhat)^2 ) )
## 
## 
## 
## 
## ###############
## 
## 
## 
## 
## ####### Let's dig in
## ####### We can use train W to construct test Xstar
## ####### Check which columns drive our response, y
## 
## X_test <- Xall[ xndx_test, ]
## y_test <- yall[ xndx_test ]
## 
## Xstar_test <- f_mactivate(U=X_test, W=xxls_out[[ length(xxls_out) ]][[ "What" ]])
## Xi <- cbind(X_test, Xstar_test)
## xlm <- lm(y_test ~ Xi)
## 
## sumxlm <- summary(xlm)
## 


###################################################
### code chunk number 3: mactivation_examples_02.Snw:340-341
###################################################
load(file="xlm4summary.RData")


###################################################
### code chunk number 4: mactivation_examples_02.Snw:347-348
###################################################
print(sumxlm)


###################################################
### code chunk number 5: b1 (eval = FALSE)
###################################################
## 
## ########### Looks like cols 1,2,4,6 are where all the action is
## ########### (in other word, mactivate didn't track signal on passes 3 and 5)
## 
## bWhat <- xxls_out[[ length(xxls_out) ]][[ "What" ]][ ,  c(1,2,4,6) ]
## bWhat
## 
## wwmag <- apply(bWhat, 1, function(x) { return(sum(abs(x)))} ) ; wwmag
## 
## plot(wwmag, type="h", lwd=4,
## main="W Coefficient Total Magnitute vs Input Term",
## xlab="Column of U (in this case, same as X)",
## ylab="Sum of magnitudes in fitted W",
## cex.lab=1.3
## )
## 


###################################################
### code chunk number 6: c1 (eval = FALSE)
###################################################
## 
## ###########################
## 
## ########## these are the terms in X that
## ########## appear to have a polynomial contribution
## ########## to our response, y
## 
## xxhigh <- which(wwmag > 0.5) ; xxhigh
## 
## ######################
## 
## 


###################################################
### code chunk number 7: d1 (eval = FALSE)
###################################################
## 
## ###### no need to show this
## 
## ######## we can try refitting using only contributing inputs
## 
## 
## m_tot <- 8
## 
## X_train <- Xall[ xndx_train, ]
## y_train <- yall[ xndx_train ]
## 
## xcmact_hybrid <-
## f_control_mactivate(
## param_sensitivity = 10^12,
## bool_free_w = TRUE,
## w0_seed = 0.01,
## w_col_search = "alternate",
## max_internal_iter = 500,
## ss_stop = 10^(-14),
## escape_rate = 1.001,
## Wadj = 1/1,
## force_tries = 0,
## lambda = 0/10000,
## tol = 10^(-14)
## )
## 
## 
## xxnow <- Sys.time()
## xxls_out_lean <-
## f_fit_hybrid_01(
## X = X_train,
## y = y_train,
## m_tot = m_tot,
## U = X_train[ , xxhigh ], ### only consider 'significant' contributors
## m_start = 1,
## mact_control = xcmact_hybrid,
## verbosity = 1
## )
## cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )
## 
## 
## 


###################################################
### code chunk number 8: aa1 (eval = FALSE)
###################################################
## 
## 
## data(df_hospitals_ortho)
## 
## xvars <- attr(df_hospitals_ortho, "modelvars")
## 
## xmx <- as.matrix(df_hospitals_ortho[ , xvars])
## 
## old_par <- par(mfrow=c(3,4))
## 
## for(jj in 1:ncol(xmx)) {
##     hist(xmx[ ,jj], main=colnames(xmx)[jj])
## }
## ######## pause and look
## 
## ################# transform
## 
## xlog_list <-
## c(
## "tot_sales",
## "tot_knee",
## "tot_hip",
## "beds",
## "rehab_beds",
## "outpatient_visits",
## "adm_costs",
## "revenue_inpatient"
## )
## 
## 
## ymx <- xmx
## 
## for(jj in 1:ncol(xmx)) {
##     ymx[ ,jj] <- log( xmx[ ,jj] + 1 )
## }
## 
## 
## for(jj in 1:ncol(ymx)) {
##     hist(ymx[ ,jj], main=colnames(ymx)[jj])
## }
## ######## pause and look
## 
## 
## 
## #### standardize
## 
## ymx_stnd <- t( ( t( ymx ) - apply(ymx, 2, mean)) / apply(ymx, 2, sd) )
## 
## 
## for(jj in 1:ncol(ymx_stnd)) {
##     hist(ymx_stnd[ ,jj], main=colnames(ymx_stnd)[jj])
## }
## ######## pause and look
## 
## 
## 
## ############# let's fit an MLR model (in place, no train-test)
## ydf_stnd <- as.data.frame(ymx_stnd)
## 
## xlm <- lm( tot_sales ~ . , data=ydf_stnd )
## 
## summary(xlm)
## 
## yhat <- xlm$fitted
## yy <- ydf_stnd[ , "tot_sales" ]
## 
## rmse <- sqrt( mean( (yy - yhat)^2 ) ) ; rmse
## 
## 1 - rmse^2 / 1  #### r-squared
## 
## 
## ######## now let's break out m-activation using 10-fold CV
## 
## 
## library(mactivate)
## 
## yall <- ymx_stnd[ , "tot_sales" ]
## Xall <- ymx_stnd[ , -1 ]
## Uall <- Xall
## 
## xfolds <- rep_len( 1:10, length(yall) ) ## 10-fold CV
## 
## m_tot <- 5
## 
## xmcont <-
## f_control_mactivate(
## param_sensitivity = 10^11,
## w0_seed = 0.1,
## max_internal_iter = 500,
## w_col_search = "one",
## bool_headStart = FALSE,
## ss_stop = 10^(-11),
## escape_rate = 1.001,
## step_size = 1/100,
## Wadj = 1/1,
## force_tries = 0,
## lambda = 0,
## tol = 10^(-8)
## )
## 
## 
## ufolds <- sort(unique(xfolds))
## 
## yout <- numeric(length(yall))
## 
## #### takes about 5-10 minutes
## 
## for(iif in 1:length(ufolds)) {
##     
##     xmask_fold <- xfolds %in% ufolds[ iif ]
##     
##     xxls_out <-
##     f_fit_hybrid_01(
##     X=Xall[ !xmask_fold, ],
##     y=yall[ !xmask_fold ],
##     m_tot=m_tot,
##     U = Xall[ !xmask_fold, ],
##     m_start = 1,
##     mact_control = xmcont,
##     verbosity = 0
##     )
##     
##     xxls_out
##     
##     yhatG <- predict(
## 	object=xxls_out, 
## 	X0=Xall[ xmask_fold, ], 
## 	U0=Uall[ xmask_fold, ],
## 	mcols=m_tot 
## 	)
##     
##     yout[ xmask_fold ] <- yhatG
##     
##     cat("Done this fold:", iif, "\n\n")
##     
## }
## 
## mact_rmse <- sqrt( mean( (yall  -  yout)^2 ) ) ; mact_rmse
## 
## cat("TT R2:", 1 - mact_rmse^2 / 1, "\n") #### m-activation R^2
## 
## par(old_par)
## 
## 


###################################################
### code chunk number 9: bb1 (eval = FALSE)
###################################################
## 
## 
## library(mactivate)
## 
## set.seed(777)
## 
## d <- 25
## N <- 100000
## 
## X <- matrix(rnorm(N*d, 0, 1), N, d) ####
## 
## colnames(X) <- paste0("x", I(1:d))
## 
## ############# primary effects
## b <- rep_len( c(-1/4, 1/4), d )
## 
## 
## ystar <-
## X %*% b +
## 1/3 * (X[ , 1]+1) * (X[ , 2]-1) * (X[ , 3]+1) -
## 1/2 * (X[ , 3]+0) * (X[ , 4]+1) * (X[ , 5]-0) * (X[ , 6]-1) +
## 1/3 * (X[ , 6]+1) * (X[ , 7]-1) -
## 1/2 * (X[ , 1]-1) * (X[ , 2]+0) * (X[ , 3]+1) * (X[ , 4]-0) * 
## (X[ , 5]-0) * (X[ , 7]-1)
## 
## m_tot <- 10
## #############
## 
## 
## xs1 <- "y ~ . + x1:x2:x3 + x3:x4:x5:x6 + x6:x7 + x1:x2:x3:x4:x5:x7"
## xs2 <- "y ~ . + x1*x2*x3 + x3*x4*x5*x6 + x6*x7 + x1*x2*x3*x4*x5*x7"
## 
## xnotQuiteTrue_formula <- eval(parse(text=xs1))
## xtrue_formula <- eval(parse(text=xs2))
## xnoint_formula <- eval(parse(text="y ~ ."))
## 
## 
## 
## 
## ysigmoid <- 1 / (1 + exp(-ystar))
## 
## range(ysigmoid)
## 
## y <- rbinom(size=1, n=N, prob=ysigmoid)
## 
## ########## standardize X
## Xall <- t( ( t(X) - apply(X, 2, mean) ) / apply(X, 2, sd) )
## yall <- y
## Nall <- N
## 
## 
## ####### fold index
## xxfoldNumber <- rep_len(1:2, N)
## 
## ufolds <- sort(unique(xxfoldNumber)) ; ufolds
## 
## 
## ############### predict
## ############### predict
## 
## dfx <- data.frame("y"=yall, Xall)
## 
## ###################
## 
## xglm <- glm(xnoint_formula , data=dfx, family=binomial(link="logit"))
## summary(xglm)
## yhat <- predict(xglm, newdata=dfx, type="response")
## mean( f_logit_cost(y=yall, yhat=yhat) )
## 
## 
## ####### known true when zero centered
## xglm <- glm(xnotQuiteTrue_formula , data=dfx, family=binomial(link="logit"))
## summary(xglm)
## yhat <- predict(xglm, newdata=dfx, type="response")
## mean( f_logit_cost(y=yall, yhat=yhat) )
## 
## 
## ####### known true when not zero centered
## xglm <- glm(xtrue_formula , data=dfx, family=binomial(link="logit"))
## summary(xglm)
## yhat <- predict(xglm, newdata=dfx, type="response")
## mean( f_logit_cost(y=yall, yhat=yhat) )
## 
## 
## 
## ################################ alternate
## ################################ alternate -- about 1.5 hours
## 
## xcmact_gradient <-
## f_control_mactivate(
## param_sensitivity = 10^9,
## bool_free_w = TRUE,
## w0_seed = 0.05,
## w_col_search = "alternate",  
## bool_headStart = TRUE,
## ss_stop = 10^(-9), ###
## escape_rate = 1.003,
## Wadj = 1/1,
## force_tries = 0
## )
## 
## 
## Uall <- Xall
## 
## xthis_fold <- ufolds[ 1 ]
## 
## xndx_test <- which( xxfoldNumber %in% xthis_fold )
## xndx_train <- setdiff( 1:Nall, xndx_test )
## 
## X_train <- Xall[ xndx_train, , drop=FALSE ]
## y_train <- yall[ xndx_train ]
## 
## xxnow <- Sys.time()
## xxls_out <-
## f_fit_gradient_logistic_01(
## X = X_train,
## y = y_train,
## m_tot = m_tot,
## U = X_train,
## m_start = 1,
## mact_control = xcmact_gradient,
## verbosity = 1
## )
## cat( difftime(Sys.time(), xxnow, units="mins"), "\n" )
## 
## 
## 
## U_test <- Xall[ xndx_test, , drop=FALSE ]
## X_test <- Xall[ xndx_test, , drop=FALSE ]
## y_test <- yall[ xndx_test ]
## 
## 
## yhatTT <- matrix(NA, length(xndx_test), m_tot+1)
## 
## for(iimm in 0:m_tot) {
##     yhat_fold <- predict(object=xxls_out, X0=X_test, U0=U_test, mcols=iimm )
##     yhatTT[ , iimm + 1 ] <- yhat_fold[[ "p0hat" ]]
## }
## 
## errs_by_m <- NULL
## for(iimm in 1:ncol(yhatTT)) {
##     yhatX <- yhatTT[ , iimm]
##     errs_by_m[ iimm ] <- mean( f_logit_cost(y=y_test, yhat=yhatX) )
##     cat(iimm, "::", errs_by_m[ iimm ])
## }
## 
## 
## ##### plot test Logit vs m
## 
## plot(0:(length(errs_by_m)-1), errs_by_m, type="l", xlab="m", ylab="Logit Cost")
## 
## 
## 
## ################## use known 'true' in glm()
## 
## xglm <- glm(xtrue_formula , data=dfx[ xndx_train, ], family=binomial(link="logit"))
## yhat <- predict(xglm, newdata=dfx[ xndx_test, ], type="response")
## 
## mean( f_logit_cost(y=y_test, yhat=yhat) )
## 
## 


