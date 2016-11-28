#' Validate and preprocess the data
#'
#' @param      dat   The dat
#' @param      dW    The d w
#' @param      nbin  how many levels should continuous W's be binned into
#'
#' @return
#' @export
#'
check_and_preprocess <- function(dat, dW, nbin = 4) {
    to_keep <- (dat$T.tilde != 0) & (dat$T.tilde != max(dat$T.tilde))
    dW <- dW[to_keep]
    dat <- dat[to_keep,]

    n.data <- nrow(dat)
    W_names <- grep('W', colnames(dat), value = TRUE)
    # ==================================================================================
    # input validation
    # ==================================================================================
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
    # ==================================================================================
    # continuous W: binning
    # ==================================================================================
    checkBinary <- function(v, naVal="NA") {
        if (!is.numeric(v)) stop("Only numeric vectors are accepted.")

        vSet = unique(v)
        if (!missing(naVal)) vSet[vSet == naVal] = NA
        vSet = vSet[!is.na(vSet)]

        if (any(as.integer(vSet) != vSet)) "con"
        else if (length(vSet) > 2) "con"
        else "bin"
    }
    checkBinary_df <- function(df, W_names) {
        all_out <- character()
        for (it in W_names) {
            all_out <- c(all_out, checkBinary(v = df[,it]))
        }
        return(all_out)
    }
    is_conti <- checkBinary_df(df = dat, W_names = W_names) != 'bin'
    W_conti <- W_names[is_conti]

    if(any(is_conti)) message(paste("Binning continuous covariates:", W_conti, collapse = ','))
    for (it in W_conti) {
        dat[,it] <- as.numeric(cut(dat[,it], nbin))
    }
    # ==================================================================================
    # check for positivity
    # ==================================================================================
    check_positivity(dat = dat, posit_level = 0.05)

    return(list(dat = dat, dW = dW, n.data = n.data, W_names = W_names))
}