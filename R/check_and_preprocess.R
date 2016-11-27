##
## @brief      validate and preprocess the data
##
## @param      dat   The dat
## @param      dW    The d w
##
## @return
##
check_and_preprocess <- function(dat, dW) {

    to_keep <- (dat$T.tilde != 0) & (dat$T.tilde != max(dat$T.tilde))
    dW <- dW[to_keep]
    dat <- dat[to_keep,]

    n.data <- nrow(dat)

    W_names <- grep('W', colnames(dat), value = TRUE)
    # ================================================================================================
    # input validation
    # ================================================================================================
    if (length(dW) != n.data) {
        stop('The length of input dW is not same as the sample size!')
    }

    if (!('delta' %in% colnames(dat))) {
        warning('delta not found. Set delta = 1.')
        dat$delta <- rep(1, nrow(dat))
    }

    if ('T.TILDE' %in% toupper(colnames(dat))) {
        # if there is t.tilde in variable, remove any T
        keeps <- setdiff(colnames(dat), 'T')
        dat <- dat[,keeps]
    }else if('T' %in% toupper(colnames(dat))){
        message('no t.tilde, rename T to T.tilde')
        dat <- rename(dat, T.tilde = T)
    }else{
        # if there are no T.tilde
        stop("There should be T.tilde variable!")
    }


    return(list(dat = dat, dW = dW, n.data = n.data, W_names = W_names))
}