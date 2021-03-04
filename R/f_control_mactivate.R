f_control_mactivate <-
function (param_sensitivity = 10^9, bool_free_w = FALSE, w0_seed = 0.10000000000000001, 
    max_internal_iter = 500, w_col_search = "one", bool_headStart = FALSE, 
    antifreeze = FALSE, ss_stop = 10^(-8), escape_rate = 1.004, 
    step_size = 1/100, Wadj = 1/1, force_tries = 0, lambda = 0, 
    tol = 10^(-8)) 
{
    if (w_col_search == "one") {
        bool_fix_w <- TRUE
        bool_alt_w <- FALSE
    }
    if (w_col_search == "all") {
        bool_fix_w <- FALSE
        bool_alt_w <- FALSE
    }
    if (w_col_search == "alternate") {
        bool_fix_w <- TRUE
        bool_alt_w <- TRUE
    }
    xls_out <- list()
    xls_out[["param_sensitivity"]] <- param_sensitivity
    xls_out[["bool_free_w"]] <- bool_free_w
    xls_out[["w0_seed"]] <- w0_seed
    xls_out[["bool_fix_w"]] <- bool_fix_w
    xls_out[["bool_alt_w"]] <- bool_alt_w
    xls_out[["bool_headStart"]] <- bool_headStart
    xls_out[["antifreeze"]] <- antifreeze
    xls_out[["ss_stop"]] <- ss_stop
    xls_out[["escape_rate"]] <- escape_rate
    xls_out[["step_size"]] <- step_size
    xls_out[["Wadj"]] <- Wadj
    xls_out[["force_tries"]] <- force_tries
    xls_out[["tol"]] <- tol
    xls_out[["max_internal_iter"]] <- max_internal_iter
    xls_out[["lambda"]] <- lambda
    class(xls_out) <- c(class(xls_out), "control_mactivate_obj")
    return(xls_out)
}
