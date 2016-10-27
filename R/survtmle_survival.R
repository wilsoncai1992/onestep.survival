#' Title
#'
#' @param dat
#'
#' @return
#' @export
#'
#' @examples
#' @import survtmle
#' @import Matrix
survtmle_survival <- function(dat) {
    T.uniq <- unique(sort(dat$T))
    # remove the first and last time points when algorithm fail
    T.uniq <- T.uniq[-1]
    T.uniq <- T.uniq[-length(T.uniq)]

    # create function inputs
    # remove failure time = 0
    ftime <- dat$T[dat$T!=0]
    if ('delta' %in% colnames(dat)) {
        ftype <- dat$delta[dat$T!=0]
    }else{
        # no censoring in the dataset
        # censoring to be all 1
        ftype <- rep(1, length(ftime))
    }
    trt <- dat$A[dat$T!=0]
    # get all W_ covariates
    W_name <- grep(names(dat), pattern = 'W', value = TRUE)
    adjustVars <- as.data.frame(dat[dat$T!=0,W_name])


    # ====================================================================================
    # compute values for all time points
    # ====================================================================================
    fit_max_time <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                             t0 = max(T.uniq),
                             SL.ftime = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                             SL.ctime = c("SL.glm","SL.mean"),
                             SL.trt = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                             # glm.ftime = paste(c('trt', W_name), collapse = ' + '),
                             # glm.trt = paste(W_name, collapse = ' + '),
                             method="hazard", returnModels = TRUE)
    # 7.8min
    allTimes <- timepoints(object = fit_max_time, times = T.uniq, returnModels = FALSE)

    s_vec <- sapply(allTimes, function(x) 1 - x$est['1 1',])
    s_df <- data.frame(s_vec, T.uniq)
    return(survival_df = s_df)
}
