#' wrapper for censoring regression SuperLearner
#'
#' @param dat data.frame, with col: T, A, W, Delta
#' @param T.uniq the unique time points of failure
#' @param Delta.SL.Lib library for SuperLearner
#'
#' @return h.hat.t matrix, nrow = n.sample, ncol = length(T.uniq);
#' each row is the predicted conditional hazard h(T=t | T>=t, A, W) for that subject
#' @export
#'
#' @examples
#' # TO DO
#' @import SuperLearner
#' @import survtmle
#' @import dplyr
#' @import tidyr
censor_SL_wrapper <- function(dat,
                              T.uniq,
                              Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")
                              ) {
    # transform original data into SL-friendly format
    dat_david <- dat

    dat_david <- rename(dat_david, ftime = T.tilde)
    dat_david <- rename(dat_david, trt = A)

    if ('ID' %in% toupper(colnames(dat_david))) {
        # if there are already id in the dataset
        dat_david <- rename(dat_david, id = ID)
    }else{
        # if no id exist
        # create 'id' on our own
        dat_david$id <- 1:nrow(dat_david)
    }

    # censoring
    if ('delta' %in% colnames(dat_david)) {
        dat_david <- rename(dat_david, ftype = delta)
    }

    # remove all other useless columns
    baseline_name <- grep('W', colnames(dat_david), value = TRUE)
    keeps <- c("id", baseline_name, 'ftime', 'ftype', 'trt')
    dat_david <- dat_david[, keeps]

    T.uniq <- unique(sort(dat_david$ftime))
    T.max <- max(T.uniq)

    adjustVars <- dat_david[,baseline_name]
    # ====================================================================================================
    # dataList
    # ====================================================================================================
    datalist <- survtmle:::makeDataList(dat = dat_david,
                                       J = 1, # one kind of failure
                                       ntrt = 2, # one kind of treatment
                                       uniqtrt = c(0,1),
                                       t0 = T.max, # time to predict on
                                       bounds=NULL)
    # ====================================================================================================
    # estimate g_2 (censoring)
    # perform censoring SL for the maximum time point
    # ====================================================================================================
    g2_hat <- survtmle:::estimateCensoring(dataList = datalist, adjustVars = adjustVars,
                                t0 = T.max,
                                ntrt = 2, # one kind of treatment
                                uniqtrt = c(0,1),
                                SL.ctime = Delta.SL.Lib,
                                returnModels = FALSE,verbose = FALSE)
    dataList2 <- g2_hat$dataList
    # ----------------------------------------------------------------------------------------------------
    # turn into wide format
    out_censor <- dataList2$`1`
    out_censor <- out_censor[,c('id', 't', 'G_dC')]
    out_censor <- spread(out_censor, t, G_dC)
    # the colname number correspond to h_{T>t-1 | T>=t-1}
    rownames(out_censor) <- out_censor$id

    # turn NA entries (after failure) into zero hazard
    out_censor[is.na(out_censor)] <- 0

    # remove the id column
    out_censor_2 <- out_censor[,-1]

    # subset the columns for those only in T.uniq
    out_censor <- out_censor_2[,T.uniq]

    return(list(out_censor = out_censor,
                out_censor_full = out_censor_2))
}
