#' Perform one-step TMLE update of survival curve
#'
#' @param D1.t.func.prev n*p matrix of previous influence curve
#' @param Pn.D1.func.prev p vector of previous mean influence curve
#' @param dat input data.frame
#' @param T.uniq grid of unique event times
#' @param W_names vector of the names of baseline covariates
#' @param dW dynamic intervention
#'
#' @return
#' @export
#'
#' @examples
#' # TO DO
#' @importFrom dplyr left_join
compute.update <- function(D1.t.func.prev, Pn.D1.func.prev, dat, T.uniq, W_names, dW) {
    # formula on p.30
    # result <- l2.inner.step(Pn.D1.t, D1.t, T.uniq) /
    # sqrt(l2.inner.step(D1.t, D1.t, T.uniq))
    # formula on p.28
    # result <- l2.inner.step(Pn.D1.func.prev, D1.t.func.prev, T.uniq) /
    # sqrt(l2.inner.step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))

    # WILSON MADE: MAY BE WRONG
    # result <- l2.inner.step(abs(Pn.D1.func.prev), D1.t.func.prev, T.uniq) /
    # sqrt(l2.inner.step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))

    # WILSON 2: MAY BE WRONG
    # calculate the number inside exp{} expression in submodel
    # numerator <- sweep(D1.t.func.prev, MARGIN=2, abs(Pn.D1.func.prev),`*`)
    # result <- numerator /
    # sqrt(l2.inner.step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))

    # ORIGINAL PAPER
    # calculate the number inside exp{} expression in submodel
    # each strata of Q is updated the same

    # numerator <- sweep(D1.t.func.prev, MARGIN=2, -abs(Pn.D1.func.prev),`*`)
    numerator <- sweep(D1.t.func.prev, MARGIN=2, Pn.D1.func.prev,`*`)
    result <- numerator /
        sqrt(l2.inner.step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))

    strata <- data.frame(A = dat$A, W = dat[,W_names])
    names(strata) <- c('A', W_names)
    colnames(result) <- paste('X', 1:ncol(result), sep = '')
    result2 <- cbind(strata, result)

    Xname <- paste(colnames(result), collapse = ' ,')
    # calculate the update value for each strata
    # eval the following command automatically with all columns in the result matrix
    # df2=aggregate(cbind(x1, x2)~A+W, data=result2, sum, na.rm=TRUE)
    eval(parse(text = paste('df2=aggregate(cbind(' , Xname, ')~A+',paste(W_names, collapse = ' +') , ', data=result2, sum, na.rm=TRUE)')))

    # MAY FAIL
    # 09-18: also update those who are not A == dW
    strata[dat$A != dW,'A'] <- dW[dat$A != dW]

    # assign the update value to each unique strata. within each strata, all update value are the same
    result_new <- left_join(strata, df2, by=c("A",W_names))
    # remove the strata dummies
    eval(
        parse(text = paste('result_new <- as.matrix(subset(result_new, select=-c(A,', paste(W_names, collapse = ','), ')))'))
    )

    # 2016-10-05: adjust back to ORIGINAL PAPER method
    one_col <- compute.step.cdf(pdf.mat = result_new, t.vec = T.uniq, start = Inf)[,1]
    result_new <- matrix(one_col,nrow = nrow(dat),ncol = ncol(result_new))

    # if (is.na(result)) {
    # result <- NA
    # }
    # result <- as.vector(result)
    # return(result)
    return(result_new)
}
