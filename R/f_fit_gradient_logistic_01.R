f_fit_gradient_logistic_01 <-
function (X, y, m_tot, U = NULL, m_start = 1, mact_control = f_control_mactivate(), 
    verbosity = 2) 
{
    w0_seed <- mact_control[["w0_seed"]]
    xbool_free_w <- mact_control[["bool_free_w"]]
    xparam_sensitivity <- mact_control[["param_sensitivity"]]
    xbool_fix_w <- mact_control[["bool_fix_w"]]
    xbool_alt_w <- mact_control[["bool_alt_w"]]
    xbool_headStart <- mact_control[["bool_headStart"]]
    xss_stop <- mact_control[["ss_stop"]]
    xstep_size <- mact_control[["step_size"]]
    xescape_rate <- mact_control[["escape_rate"]]
    xWadj <- mact_control[["Wadj"]]
    xforce_tries <- mact_control[["force_tries"]]
    if (is.null(U)) {
        U <- X
    }
    d <- ncol(X)
    N <- nrow(X)
    du <- ncol(U)
    xls_out <- list()
    Xint <- cbind(int = rep(1, N), X)
    xglm <- glm(y ~ Xint - 1, family = binomial(link = "logit"))
    bhatX <- xglm[["coefficients"]]
    bhatX
    bhatX <- unname(bhatX)
    bhatX
    bhatX[is.na(bhatX)] <- 0
    xls_out[[1]] <- list(What = matrix(0, d, 0), cchat = numeric(0), 
        b0hat = bhatX[1], bbhat = bhatX[2:(d + 1)])
    if (xbool_headStart) {
        ib0 <- bhatX[1]
        ibb <- bhatX[2:(d + 1)]
    }
    else {
        ib0 <- 0
        ibb <- rep(0, d)
    }
    icc <- rep(0, 1)
    iW <- matrix(w0_seed, du, 1)
    iXstar <- matrix(0, N, 1)
    iim <- 1
    m_start <- 1
    for (iim in m_start:m_tot) {
        if (iim > 1) {
            iW0 <- matrix(w0_seed, du, iim)
            iW0[, 1:(iim - 1)] <- iW[, 1:(iim - 1)]
            iW <- iW0
            iXstar0 <- matrix(0, N, iim)
            iXstar0[, 1:(iim - 1)] <- iXstar[, 1:(iim - 1)]
            iXstar <- iXstar0
            icc0 <- rep(0, iim)
            icc0[1:(iim - 1)] <- icc[1:(iim - 1)]
            icc <- icc0
        }
        rownames(iW) <- make.names(colnames(U), unique = TRUE, 
            allow_ = TRUE)
        if (iim > 1 & xbool_alt_w) {
            jjj_rng <- c(1, 2)
        }
        else {
            jjj_rng <- c(1)
        }
        for (jjj in jjj_rng) {
            if (jjj == 1) {
                xbool_fix_w_use <- xbool_fix_w
            }
            else {
                xbool_fix_w_use <- FALSE
            }
            xdeltaCO <- 1/xparam_sensitivity
            xstep_size_use <- xstep_size
            xbool_keep_going <- TRUE
            kk_i <- 0
            while (xbool_keep_going) {
                xdeltaCO <- xdeltaCO * xescape_rate
                if (xbool_fix_w_use) {
                  iXstar[, iim] <- f_mactivate(U = U, W = iW[, 
                    iim, drop = FALSE])
                }
                else {
                  iXstar <- f_mactivate(U = U, W = iW)
                }
                iyhat <- ib0 + X %*% ibb + iXstar %*% icc
                e_niyhat <- exp(-iyhat)
                df_dyhat <- -y * e_niyhat/(e_niyhat + 1) + (y - 
                  1) * e_niyhat/((1/(e_niyhat + 1) - 1) * (e_niyhat + 
                  1)^2)
                df_dyhat[is.na(df_dyhat) | is.nan(df_dyhat)] <- 0
                iphat <- 1/(e_niyhat + 1)
                ierrs <- numeric(length(iphat))
                ierrs[y == 0] <- -2 * log(1 - iphat[y == 0])
                ierrs[y == 1] <- -2 * log(iphat[y == 1])
                df_db0 <- as.vector(2 * sum(df_dyhat)/N)
                df_db0
                df_dbb <- as.vector(2 * t(df_dyhat) %*% X/N)
                df_dbb
                df_dcc <- as.vector(2 * t(df_dyhat) %*% iXstar/N)
                df_dcc
                if (xbool_fix_w_use) {
                  df_dW <- f_dmss_dW(U = U, Xstar = iXstar[, 
                    iim, drop = FALSE], W = iW[, iim, drop = FALSE], 
                    yerrs = df_dyhat, cc = icc[iim])/N
                }
                else {
                  df_dW <- f_dmss_dW(U = U, Xstar = iXstar, W = iW, 
                    yerrs = df_dyhat, cc = icc)/N
                }
                ib0p <- ib0
                ibbp <- ibb
                iccp <- icc
                iWp <- iW
                xbool_stepTry <- TRUE
                while (xbool_stepTry & xstep_size_use > xss_stop) {
                  ib0 <- ib0 - xstep_size_use * df_db0
                  ibb <- ibb - xstep_size_use * df_dbb
                  if (xbool_fix_w_use) {
                    iW[, iim] <- iW[, iim] - xstep_size_use * 
                      df_dW * xWadj
                    icc[iim] <- icc[iim] - xstep_size_use * df_dcc[iim]
                  }
                  else {
                    iW <- iW - xstep_size_use * df_dW * xWadj
                    icc <- icc - xstep_size_use * df_dcc
                  }
                  if (!xbool_free_w) {
                    iW[iW < 0] <- 0
                    iW[iW > 1] <- 1
                  }
                  if (xbool_fix_w_use) {
                    iXstar[, iim] <- f_mactivate(U = U, W = iW[, 
                      iim, drop = FALSE])
                  }
                  else {
                    iXstar <- f_mactivate(U = U, W = iW)
                  }
                  iyhat <- ib0 + X %*% ibb + iXstar %*% icc
                  iphat <- 1/(exp(-iyhat) + 1)
                  post_ierrs <- numeric(length(iphat))
                  post_ierrs[y == 0] <- -2 * log(1 - iphat[y == 
                    0])
                  post_ierrs[y == 1] <- -2 * log(iphat[y == 1])
                  xpre_rmse <- (mean(ierrs))
                  xpre_rmse
                  xpost_rmse <- (mean(post_ierrs))
                  xpost_rmse
                  xpre_rmse - xpost_rmse
                  if (10^6 * (xpre_rmse - xpost_rmse) < xpre_rmse * 
                    xdeltaCO) {
                    ib0 <- ib0p
                    ibb <- ibbp
                    icc <- iccp
                    iW <- iWp
                    xstep_size_use <- xstep_size_use/3
                  }
                  else {
                    xstep_size_use <- xstep_size_use * 1.3
                    xbool_stepTry <- FALSE
                  }
                }
                if (xstep_size_use <= xss_stop & kk_i >= xforce_tries) {
                  xbool_keep_going <- FALSE
                }
                kk_i <- kk_i + 1
                xxcatOut <- ""
                if (verbosity >= 3) {
                  xxcatOut <- paste0("-- log CO: ", log(xdeltaCO, 
                    10))
                }
                if (verbosity >= 1) {
                  cat(kk_i, xpre_rmse, xpost_rmse, "--", xstep_size_use, 
                    "-- pre m post delta:", xpre_rmse - xpost_rmse, 
                    xxcatOut, "--", "cc:", icc, "\n")
                }
            }
            cat("found m =", iim, " -- ", jjj, "\n")
            print(iW)
            print(icc)
        }
        xls_out[[iim + 1]] <- list(What = iW, cchat = icc, b0hat = ib0, 
            bbhat = ibb)
    }
    class(xls_out) <- c(class(xls_out), "mactivate_fit_gradient_logistic_01")
    return(xls_out)
}
