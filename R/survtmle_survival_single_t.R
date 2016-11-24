#' Iterative TMLE for survival curve at specific time point
#'
#' @param dat data.frame with columns T, A, C, W. All columns with character "W" will be treated as baseline covariates.
#' @param tk time point to compute survival probability
#'
#' @return
#' @export
#'
#' @examples
#' @import survtmle
#' @import Matrix
survtmle_survival_single_t <- function(dat, tk,
                                       dW = rep(1, nrow(dat))
                                       ) {
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
        # no censoring in the dataset
        # censoring to be all 1
        ftype <- rep(1, length(ftime))
    }
    # since survtmle can only accept counterfactual A = 1
    # we set trt = I{A = dW}, so that A' = dW
    trt <- (dat$A[dat$T.tilde!=0] == dW) + 0
    # get all W_ covariates
    W_name <- grep(names(dat), pattern = 'W', value = TRUE)
    adjustVars <- as.data.frame(dat[dat$T.tilde!=0,W_name])


    # ====================================================================================
    # compute values for all time points
    # ====================================================================================
    fit_max_time <- survtmle(ftime = ftime, ftype = ftype, trt = trt, adjustVars = adjustVars,
                             t0 = tk,
                             SL.ftime = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                             SL.ctime = c("SL.glm","SL.mean"),
                             SL.trt = c("SL.glm","SL.mean","SL.step", "SL.earth"),
                             # glm.ftime = paste(c('trt', W_name), collapse = ' + '),
                             # glm.trt = paste(W_name, collapse = ' + '),
                             method="hazard", returnModels = TRUE)
    # 7.8min
    # allTimes <- timepoints(object = fit_max_time, times = T.uniq, returnModels = FALSE)
    est <- 1 - fit_max_time$est['1 1',]
    var <- fit_max_time$var['1 1', '1 1']

    return(list(est = est,
                var = var,
                meanIC = fit_max_time$meanIC,
                ic = fit_max_time$ic))
}