#' Iterative TMLE for survival curve
#'
#' @param dat data.frame with columns T, A, W. All columns with character "W" will be treated as baseline covariates.
#'
#' @return data.frame, where the first column is survival probability, second column is the time point of the survival curve
#' @export
#'
#' @examples
#' @import survtmle
#' @import Matrix
survtmle_survival <- function(dat, dW = rep(1, nrow(dat)),
                              T.cutoff = NULL,
                              SL.ftime = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                              SL.ctime = c("SL.glm","SL.mean"),
                              SL.trt = c("SL.glm","SL.mean","SL.step", "SL.earth")) {
    # ===================================================================================
    # preparation
    # ===================================================================================
    after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
    dat <- after_check$dat
    dW <- after_check$dW
    n.data <- after_check$n.data
    W_names <- after_check$W_names

    T.uniq <- unique(sort(dat$T.tilde))

    # create function inputs
    ftime <- dat$T.tilde
    ftype <- dat$delta

    if(all(dW == 0)) {
        trt <- 1 - dat$A # when dW is all zero
    }else if(all(dW == 1)){
        trt <- dat$A
    }else{
        stop('not implemented!')
    }

    adjustVars <- as.data.frame(dat[,W_names])
    # ====================================================================================
    # compute values for all time points
    # ====================================================================================
    fit_max_time <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                             t0 = max(T.uniq),
                             SL.ftime = SL.ftime,
                             SL.ctime = SL.ctime,
                             SL.trt = SL.trt,
                             # glm.ftime = paste(c('trt', W_names), collapse = ' + '),
                             # glm.trt = paste(W_names, collapse = ' + '),
                             method="hazard", returnModels = TRUE,
                             verbose = FALSE)
    # 7.8min
    allTimes <- timepoints(object = fit_max_time, times = T.uniq, returnModels = FALSE)

    s_vec <- sapply(allTimes, function(x) 1 - x$est['1 1',])
    survival_df <- data.frame(s_vec, T.uniq)

    class(survival_df) <- 'surv_survtmle'
    return(survival_df)
}
