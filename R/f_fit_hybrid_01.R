f_fit_hybrid_01 <-
function (X, y, m_tot, U = NULL, m_start = 1, mact_control = f_control_mactivate(), 
    verbosity = 2) 
{
    w0_seed <- mact_control[["w0_seed"]]
    xbool_free_w <- mact_control[["bool_free_w"]]
    xparam_sensitivity <- mact_control[["param_sensitivity"]]
    xbool_fix_w <- mact_control[["bool_fix_w"]]
    xbool_alt_w <- mact_control[["bool_alt_w"]]
    xmax_internal_iter <- mact_control[["max_internal_iter"]]
    xss_stop <- mact_control[["ss_stop"]]
    xstep_size <- mact_control[["step_size"]]
    xescape_rate <- mact_control[["escape_rate"]]
    xWadj <- mact_control[["Wadj"]]
    xforce_tries <- mact_control[["force_tries"]]
    xtol <- mact_control[["tol"]]
    xreg <- mact_control[["lambda"]]
    if (is.null(U)) {
        U <- X
    }
    if (is.null(xreg)) {
        xreg <- 0
    }
    d <- ncol(X)
    N <- nrow(X)
    du <- ncol(U)
    Xint <- cbind(int = rep(1, nrow(X)), X)
    tXintXint <- crossprod(Xint)
    tXintXint <- tXintXint + diag(xreg * diag(tXintXint), ncol(tXintXint))
    inv_tXintXint <- solve(tXintXint)
    inv_tXintXint_tXint <- inv_tXintXint %*% t(Xint)
    xls_out <- list()
    bhatX <- as.vector(inv_tXintXint_tXint %*% y)
    bhatX <- unname(bhatX)
    bhatX[is.na(bhatX)] <- 0
    xls_out[[1]] <- list(What = matrix(0, d, 0), cchat = numeric(0), 
        bbhat = bhatX)
    icc <- rep(0, 1)
    iW <- matrix(w0_seed, du, 1)
    iim <- 1
    m_start <- 1
    for (iim in m_start:m_tot) {
        if (iim > 1) {
            iW0 <- matrix(w0_seed, du, iim)
            iW0[, 1:(iim - 1)] <- iW[, 1:(iim - 1)]
            iW <- iW0
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
            m <- iim
            xbool_keep_going <- TRUE
            xstep_size_use <- xstep_size
            kk <- 0
            while (xbool_keep_going) {
                kk <- kk + 1
                iXstar <- f_mactivate(U = U, W = iW)
                ixsicc <- iXstar %*% icc
                y_nocw <- y - ixsicc
                bhats <- inv_tXintXint_tXint %*% y_nocw
                yhatb <- Xint %*% bhats
                yfull_prior_err <- sqrt(mean((y_nocw - yhatb)^2))
                yfull_prior_err
                y_nob <- y - yhatb
                xdeltaCO <- xdeltaCO * xescape_rate^10
                xstep_size_use <- xstep_size
                iik <- 0
                while (iik < xmax_internal_iter) {
                  iik <- iik + 1
                  xdeltaCO <- xdeltaCO * xescape_rate
                  y_nob_hat <- iXstar %*% icc
                  iyw_errs <- y_nob_hat - y_nob
                  df_dcc <- as.vector(2 * t(iyw_errs) %*% iXstar/N)
                  df_dcc
                  if (xbool_fix_w_use) {
                    df_dW <- f_dmss_dW(U = U, Xstar = iXstar[, 
                      iim, drop = FALSE], W = iW[, iim, drop = FALSE], 
                      yerrs = iyw_errs, cc = icc[iim])/N
                    df_dW
                  }
                  else {
                    df_dW <- f_dmss_dW(U = U, Xstar = iXstar, 
                      W = iW, yerrs = iyw_errs, cc = icc)/N
                    df_dW
                  }
                  iccp <- icc
                  iWp <- iW
                  df_dW[is.na(df_dW) | is.nan(df_dW)] <- 0
                  xbool_stepTry <- TRUE
                  while (xbool_stepTry & xstep_size_use > xss_stop) {
                    if (xbool_fix_w_use) {
                      iW[, iim] <- iW[, iim] - xstep_size_use * 
                        df_dW * xWadj
                      icc[iim] <- icc[iim] - xstep_size_use * 
                        df_dcc[iim]
                    }
                    else {
                      iW <- iW - xstep_size_use * df_dW * xWadj
                      icc <- icc - xstep_size_use * df_dcc
                    }
                    if (!xbool_free_w) {
                      iW[iW < 0] <- 0
                      iW[iW > 1] <- 1
                    }
                    iXstar <- f_mactivate(U = U, W = iW)
                    yy_errs <- y_nob - iXstar %*% icc
                    xpre_rmse <- sqrt(mean(iyw_errs^2))
                    xpre_rmse
                    xpost_rmse <- sqrt(mean(yy_errs^2))
                    xpost_rmse
                    xpre_rmse - xpost_rmse
                    if (10^6 * (xpre_rmse - xpost_rmse) < xpre_rmse * 
                      xdeltaCO) {
                      icc <- iccp
                      iW <- iWp
                      xstep_size_use <- xstep_size_use/3
                    }
                    else {
                      xstep_size_use <- xstep_size_use * 1.3
                      xbool_stepTry <- FALSE
                    }
                  }
                  if (xstep_size_use <= xss_stop) {
                    iik <- xmax_internal_iter + 2
                  }
                  xxcatOut <- ""
                  if (verbosity >= 3) {
                    xxcatOut <- paste0("-- log CO: ", log(xdeltaCO, 
                      10))
                  }
                  if (verbosity >= 2) {
                    cat("Gradient Step:", iik, " -- ", "RMSE:", 
                      sqrt(mean((yy_errs)^2)), xxcatOut, "\n")
                  }
                }
                yy_all_hat <- Xint %*% bhats + iXstar %*% icc
                yerrs_all <- y - yy_all_hat
                yfull_post_err <- sqrt(mean((yerrs_all)^2))
                cat("Refit step:", kk, "--", "RMSE:", yfull_post_err, 
                  "-- cc:", icc, "\n")
                if (yfull_prior_err - yfull_post_err < xtol * 
                  yfull_post_err & kk > xforce_tries) {
                  xbool_keep_going <- FALSE
                }
            }
            cat("found m =", iim, " -- ", jjj, "\n")
            print(iW)
            print(icc)
        }
        xls_out[[iim + 1]] <- list(What = iW, cchat = icc, bbhat = bhats)
    }
    class(xls_out) <- c(class(xls_out), "mactivate_fit_hybrid_01")
    return(xls_out)
}
