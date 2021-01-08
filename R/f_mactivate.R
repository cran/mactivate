f_mactivate <-
function (U, W) 
{
    d <- ncol(U)
    m <- ncol(W)
    N <- nrow(U)
    Xout <- matrix(0, N, m)
    Cout <- .C("mactivate_a", as.integer(N), as.integer(d), as.integer(m), 
        as.double(as.vector(U)), as.double(as.vector(W)), Xout = as.double(as.vector(Xout)))
    Xout <- Cout[["Xout"]]
    dim(Xout) <- c(N, m)
    colnames(Xout) <- colnames(W)
    return(Xout)
}
