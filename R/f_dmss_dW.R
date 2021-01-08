f_dmss_dW <-
function (U, Xstar, W, yerrs, cc) 
{
    d <- ncol(U)
    m <- ncol(W)
    N <- nrow(U)
    Wout <- matrix(0, d, m)
    Cout <- .C("dmss_dW_b", as.integer(N), as.integer(d), as.integer(m), 
        as.double(as.vector(U)), as.double(as.vector(W)), as.double(as.vector(Xstar)), 
        as.double(yerrs), as.double(cc), Wout = as.double(as.vector(Wout)))
    Wout <- Cout[["Wout"]]
    dim(Wout) <- c(d, m)
    colnames(Wout) <- colnames(W)
    return(Wout)
}
