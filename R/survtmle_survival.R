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
                              SL.ftime = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                              SL.ctime = c("SL.glm","SL.mean"),
                              SL.trt = c("SL.glm","SL.mean","SL.step", "SL.earth")) {
    T.uniq <- unique(sort(dat$T.tilde))
    # remove the first and last time points when algorithm fail
    T.uniq <- T.uniq[-1]
    T.uniq <- T.uniq[-length(T.uniq)]

    # create function inputs
    # remove failure time = 0
    ftime <- dat$T.tilde[dat$T.tilde!=0]
    if ('delta' %in% colnames(dat)) {
        ftype <- dat$delta[dat$T.tilde!=0]
    }else{
        warning('delta not found. Set delta = 1.')
        # no censoring in the dataset
        # censoring to be all 1
        ftype <- rep(1, length(ftime))
    }
    # trt <- dat$A[dat$T.tilde!=0]

    # subset of dW
    dW <- dW[dat$T.tilde!=0]
    # WRONG: since survtmle can only accept counterfactual A = 1
    # we set trt = I{A = dW}, so that A' = dW
    # trt <- (dat$A[dat$T.tilde!=0] == dW) + 0
    if(all(dW == 0)) {
        trt <- 1 - dat$A[dat$T.tilde!=0] # when dW is all zero
    }else if(all(dW == 1)){
        trt <- dat$A[dat$T.tilde!=0]
    }else{
        stop('not implemented!')
    }

    # get all W_ covariates
    W_name <- grep(names(dat), pattern = 'W', value = TRUE)
    adjustVars <- as.data.frame(dat[dat$T.tilde!=0,W_name])
    # ====================================================================================
    # compute values for all time points
    # ====================================================================================
    fit_max_time <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                             t0 = max(T.uniq),
                             SL.ftime = SL.ftime,
                             SL.ctime = SL.ctime,
                             SL.trt = SL.trt,
                             # glm.ftime = paste(c('trt', W_name), collapse = ' + '),
                             # glm.trt = paste(W_name, collapse = ' + '),
                             method="hazard", returnModels = TRUE)
    # 7.8min
    allTimes <- timepoints(object = fit_max_time, times = T.uniq, returnModels = FALSE)

    s_vec <- sapply(allTimes, function(x) 1 - x$est['1 1',])
    survival_df <- data.frame(s_vec, T.uniq)

    class(survival_df) <- 'survtmle_survival'
    return(survival_df)
}
