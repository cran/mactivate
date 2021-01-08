predict.mactivate_fit_hybrid_01 <-
function (object, X0, U0 = NULL, mcols, ...) 
{
    if (is.null(U0)) {
        U0 <- X0
    }
    xmcol <- mcols
    mobj_use <- object[[xmcol + 1]]
    What_use <- mobj_use[["What"]]
    X0star <- f_mactivate(U = U0, W = What_use)
    Xfull <- cbind(rep(1, nrow(X0)), X0, X0star)
    bc_hat <- c(mobj_use[["bbhat"]], mobj_use[["cchat"]])
    y0hat <- Xfull %*% bc_hat
    return(y0hat)
}
