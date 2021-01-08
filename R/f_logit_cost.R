f_logit_cost <-
function (y, yhat) 
{
    yonemask <- y == 1
    ylerrs <- numeric(length(y))
    ylerrs[yonemask] <- -2 * log(yhat[yonemask])
    ylerrs[!yonemask] <- -2 * log(1 - yhat[!yonemask])
    return(ylerrs)
}
