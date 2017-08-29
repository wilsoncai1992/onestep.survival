#' One-step TMLE estimator for survival at specific time point; Loop over all times
#'
#' @param dat data.frame with columns T, A, C, W. All columns with character "W" will be treated as baseline covariates.
#' @param tk time point to compute survival probability
#' @param dW binary input vector specifying dynamic treatment (as a function output of W)
#' @param SL.trt SuperLearner library for fitting treatment regression
#' @param SL.ctime SuperLearner library for fitting censoring regression
#' @param SL.ftime SuperLearner library for fitting conditional hazard regression
#' @param maxIter maximal number of recursion for one-step
#' @param epsilon_step step size for one-step recursion
#' @param tol tolerance for optimization
#' @param T.cutoff  manual right censor the data; remove parts dont want to esimate
#' @param verbose to print log-likelihood value during optimzation
#'
#' @return
#' @export
#'
#' @examples
#' # TO DO
#' @import dplyr
#' @import survtmle
onestep_single_all_t <- function(dat, dW = rep(1, nrow(dat)),
                             SL.trt = c("SL.glm", "SL.step", "SL.earth"),
                             SL.ctime = c("SL.glm", "SL.step", "SL.earth"),
                             SL.ftime = c("SL.glm", "SL.step", "SL.earth"),
                             maxIter = 3e2,
                             epsilon_step = 1e-3,
                             tol = 1/nrow(dat),
                             T.cutoff = NULL,
                             verbose = FALSE){
    # ====================================================================================================
    # input validation
    # ====================================================================================================
    after_check <- check_and_preprocess(dat = dat, dW = dW, T.cutoff = T.cutoff)
    dat <- after_check$dat
    dW <- after_check$dW
    n.data <- after_check$n.data
    W_names <- after_check$W_names
    # ====================================================================================================
    # preparation: make data in survtmle format (dat_david)
    # ====================================================================================================
    # transform original data into SL-friendly format
    dat_david <- dat

    dat_david <- rename(dat_david, ftime = T.tilde)
    dat_david <- rename(dat_david, trt = A)

    if(all(dW == 0)) {
        dat_david$trt <- 1 - dat_david$trt # when dW is all zero
    }else if(all(dW == 1)){

    }else{
        stop('not implemented!')
    }

    if ('ID' %in% toupper(colnames(dat_david))) {
        # if there are already id in the dataset
        dat_david <- rename(dat_david, id = ID)
    }else{
        warning('no id exist, create \'id\' on our own')
        dat_david$id <- 1:nrow(dat_david)
    }

    # censoring
    dat_david <- rename(dat_david, ftype = delta)

    # remove all other useless columns
    baseline_name <- W_names
    keeps <- c("id", baseline_name, 'ftime', 'ftype', 'trt')
    dat_david <- dat_david[,keeps]
    # ====================================================================================================
    # prepare
    # ====================================================================================================
    T.uniq <- unique(sort(dat_david$ftime))
    T.max <- max(T.uniq)

    adjustVars <- dat_david[,baseline_name, drop = FALSE]
    # ====================================================================================================
    # estimate g
    # ====================================================================================================
    message('estimating g_1')
    g1_hat <- estimateTreatment(dat = dat_david, adjustVars = adjustVars,
                                SL.trt = SL.trt,verbose = verbose, returnModels = FALSE)
    g1_dat <- g1_hat$dat
    # ====================================================================================================
    # make datalist
    # ====================================================================================================
    datalist <- survtmle::makeDataList(dat = g1_dat,
                                       J = 1, # one kind of failure
                                       ntrt = 2, # one kind of treatment
                                       uniqtrt = c(0,1),
                                       t0 = T.max, # time to predict on
                                       bounds=NULL)
    # yo <- datalist[[3]]
    # ====================================================================================================
    # estimate g_2 (censoring)
    # ====================================================================================================
    message('estimating g_2')
    g2_hat <- survtmle:::estimateCensoring(dataList = datalist, adjustVars = adjustVars,
                                t0 = T.max,
                                ntrt = 2, # one kind of treatment
                                uniqtrt = c(0,1),
                                SL.ctime = SL.ctime,
                                returnModels = FALSE,verbose = verbose)
    dataList2 <- g2_hat$dataList
    # ====================================================================================================
    # estimate h(t) (hazard)
    # ====================================================================================================
    message('estimating hazard')
    h_hat <- survtmle::estimateHazards(dataList = dataList2,
                             J = 1,
                             adjustVars = adjustVars,
                             SL.ftime = SL.ftime,
                             returnModels = FALSE,
                             verbose = verbose,
                             glm.ftime = NULL)
    dataList2 <- h_hat$dataList
    # check convergence
    suppressWarnings(if (all(dataList2[[1]] == "convergence failure")) {
        return("estimation convergence failure")
    })
    # ====================================================================================================
    # transform to survivial
    # ====================================================================================================
    dataList2 <- updateVariables(dataList = dataList2, allJ = 1,
                                 ofInterestJ = 1, nJ = 2, uniqtrt = c(0,1),
                                 ntrt = 2, t0 = T.max, verbose = verbose)
    # hehe <- dataList2[[3]]
    dataList2_before_target <- dataList2
    onestep_out_all <- list()
    tk_count <- 0
    for (tk in T.uniq) {
        tk_count <- tk_count + 1
        dataList2 <- dataList2_before_target

        # ====================================================================================================
        # get IC
        # ====================================================================================================
        dat_david2 <- getHazardInfluenceCurve(dataList = dataList2, dat = dat_david,
                                              ofInterestJ = 1, allJ = 1, nJ = 2, uniqtrt = c(0,1),
                                              ntrt = 2, verbose = verbose, t0 = tk)
        infCurves <- dat_david2[, grep("D.j", names(dat_david2))]
        meanIC <- colMeans(infCurves)
        # ====================================================================================================
        # targeting
        # ====================================================================================================
        calcLoss <- function(Y, QAW){
            -mean(Y * log(QAW) + (1-Y) * log(1 - QAW))
        }

        # if the derivative of the loss the positive, change the tergeting direction
        epsilon_step1 <- epsilon_step2 <- epsilon_step
        if (meanIC[2] < 0) { epsilon_step2 <- -epsilon_step2}
        if (meanIC[1] < 0) { epsilon_step1 <- -epsilon_step1}

        loss_old <- Inf
        # loss_new <- calcLoss(Y = dataList2$`1`$N1, QAW = dataList2$`1`$Q1Haz)
        loss_new <- calcLoss(Y = dataList2$obs$N1, QAW = dataList2$obs$Q1Haz)
        message(paste('targeting', tk))
        iter_count <- 0

        # while (any(abs(meanIC) > tol) & iter_count <= maxIter) {
        # while (any(abs(meanIC[2]) > tol) & iter_count <= maxIter) {
        while ((loss_new <= loss_old) & iter_count <= maxIter) {
            iter_count <- iter_count + 1
            if(verbose) print(loss_new)
            # print(meanIC[1,])

            # fluctuate -> update to dataList2
            dataList2$`1`$Q1Haz <- plogis(qlogis(dataList2$`1`$Q1Haz) + epsilon_step2 * dataList2$`1`$H1.jSelf.z1 + epsilon_step1 * dataList2$`1`$H1.jSelf.z0)
            dataList2$`0`$Q1Haz <- plogis(qlogis(dataList2$`0`$Q1Haz) + epsilon_step2 * dataList2$`0`$H1.jSelf.z1 + epsilon_step1 * dataList2$`0`$H1.jSelf.z0)
            dataList2$obs$Q1Haz <- plogis(qlogis(dataList2$obs$Q1Haz) + epsilon_step2 * dataList2$obs$H1.jSelf.z1 + epsilon_step1 * dataList2$obs$H1.jSelf.z0)


            # calculate survival again
            dataList2 <- updateVariables(dataList = dataList2, allJ = 1,
                                         ofInterestJ = 1, nJ = 2, uniqtrt = c(0,1),
                                         ntrt = 2, t0 = tk, verbose = verbose)
            # calculate IC again
            dat_david2 <- getHazardInfluenceCurve(dataList = dataList2, dat = dat_david2,
                                                  ofInterestJ = 1, allJ = 1, nJ = 2, uniqtrt = c(0,1),
                                                  ntrt = 2, verbose = verbose, t0 = tk)
            infCurves <- dat_david2[, grep("D.j", names(dat_david2))]
            meanIC_old <- meanIC
            meanIC <- colMeans(infCurves)

            # loss_new <- calcLoss(Y = dataList2$`1`$N1, QAW = dataList2$`1`$Q1Haz)
            loss_old <- loss_new
            loss_new <- calcLoss(Y = dataList2$obs$N1, QAW = dataList2$obs$Q1Haz)


            # if one converges, then stop update
            if ((abs(meanIC[1]) < tol) | (meanIC_old[1] * meanIC[1] <= 0)) {
                # if changes sign or converges, then stop update
                epsilon_step1 <- 0
            }
            if ((abs(meanIC[2]) < tol) | (meanIC_old[2] * meanIC[2] <= 0)) {
                # if changes sign or converges, then stop update
                epsilon_step2 <- 0
            }

            # if (all(abs(meanIC) < tol)) {
            # print(abs(meanIC))
            if (abs(meanIC)[2] < tol) {
                # all IC become zero mean
                message('Success Converge!')
                break()
            }
            if (epsilon_step1 == 0 & epsilon_step2 == 0){
                message('Success Converge! with epsilon_step too large.')
                break()
            }
        }

        if (iter_count == maxIter + 1) {
            message("TMLE fluctuations did not converge. Check that meanIC is adequately small and proceed with caution.")
        }


        # ====================================================================================================
        # get final estimates
        # ====================================================================================================
        est <- rowNames <- NULL

        # parameter estimates
        for (j in 1) {
            for (z in c(0,1)) {
                eval(parse(text = paste("est <- rbind(est, dat_david2$margF",
                                        j, ".z", z, ".t0[1])", sep = "")))
                rowNames <- c(rowNames, paste(c(z, j), collapse = " "))
            }
        }
        row.names(est) <- rowNames
        var <- t(as.matrix(infCurves)) %*% as.matrix(infCurves)/n.data^2
        row.names(var) <- colnames(var) <- rowNames

        # output static interventions
        est <- 1 - est['1 1',]
        var <- var['1 1', '1 1']

        # if (all(dW == 1)) {
        #     est <- 1 - est['1 1',]
        #     var <- var['1 1', '1 1']
        # }else if (all(dW == 0)) {
        #     est <- 1 - est['0 1',]
        #     var <- var['0 1', '0 1']
        # }
        onestep_out_all[[tk_count]] <- list(est = est, var = var, meanIC = meanIC, ic = infCurves)
    }

    s_vec <- sapply(onestep_out_all, function(x) x$est)
    survival_df <- data.frame(s_vec, T.uniq)
    class(survival_df) <- 'surv_survtmle'

    return(list(survival_df = survival_df, onestep_out_all = onestep_out_all))
}

