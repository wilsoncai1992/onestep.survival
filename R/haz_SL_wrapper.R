#' wrapper for conditonal hazard regression SuperLearner
#'
#' use T, A, C, W data format as input
#'
#' @param dat data.frame, with col: T, A, C, W
#' @param T.uniq vector of unique event times
#' @param ht.SL.Lib library for SuperLearner
#'
#' @return h.hat.t matrix, nrow = n.sample, ncol = length(T.uniq);
#'          each row is the predicted conditional hazard h(T=t | T>=t, A, W) for that subject
#' @export
#'
#' @examples
#' # TO DO
#' @importFrom plyr rename
#' @import tidyr
haz_SL_wrapper <- function(dat,
                           T.uniq,
                           ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")
                           ) {
    # ----------------------------------------------------------------------------------------------------
    # transform original data into SL-friendly format
    dat_david <- dat

    if ('T.TILDE' %in% toupper(colnames(dat_david))) {
        # if there is t.tilde in variable, remove any T
        keeps <- setdiff(colnames(dat_david), 'T.Tilde')
        dat_david <- dat_david[,keeps]
    }else if('T' %in% toupper(colnames(dat_david))){
        # if no t.tilde, rename to T.tilde
        dat_david <- plyr::rename(dat_david, c('T' = 'T.tilde'))
    }else{
        # if there are no T
        stop("There should be T variable!")
    }

    dat_david <- plyr::rename(dat_david, c('A' = 'Z'))

    if ('ID' %in% toupper(colnames(dat_david))) {
        # if there are already id in the dataset
        dat_david <- plyr::rename(dat_david, c('ID' = 'id'))
    }else{
        # if no id exist
        # create 'id' on our own
        dat_david$id <- 1:nrow(dat_david)
    }

    # censoring
    if ('delta' %in% colnames(dat_david)) {
        dat_david <- plyr::rename(dat_david, c('delta' = 'Delta.J'))
    }else{
        # no censoring in the dataset
        # censoring to be all 1
        dat_david$Delta.J <- 1
    }

    # remove all other useless columns
    baseline_name <- grep('W', colnames(dat_david), value = TRUE)
    keeps <- c("id", baseline_name, 'T.tilde', 'Delta.J', 'Z')
    dat_david <- dat_david[,keeps]

    # first fit SL on the max time point
    T.it <- max(T.uniq)
    datalist <- makeDataList(dat = dat_david,
                             J = 1, # one kind of failure
                             nZ = 2, # one kind of treatment
                             Z = c(0,1),
                             t0 = T.it, # time to predict on
                             bounds=NULL)
    # ----------------------------------------------------------------------------------------------------
    # perform hazard SL for the maximum time point
    Haz_hat <- estimateHazards(dataList = datalist, J=1, verbose=FALSE, strata=NULL,
                               adjustVars = dat_david[,baseline_name],
                               SLlibrary.event = ht.SL.Lib,
                               glmFormula.event = NULL,
                               bounds = NULL
    )
    # Haz_hat[[1]][Haz_hat[[1]]$t == T.it,'Q1Haz']
    # mean(Haz_hat[[1]][Haz_hat[[1]]$t == T.it,'Q1Haz'])

    # turn into wide format

    out_haz <- Haz_hat[[3]]
    out_haz <- out_haz[,c('id', 't', 'Q1Haz')]
    out_haz <- spread(out_haz, t, Q1Haz)
    # the colname number correspond to h_{T>t-1 | T>=t-1}
    rownames(out_haz) <- out_haz$id

    # turn NA entries (after failure) into zero hazard
    out_haz[is.na(out_haz)] <- 0

    # remove the id column
    out_haz_2 <- out_haz[,-1]

    # subset the columns for those only in T.uniq
    out_haz <- out_haz_2[,T.uniq]

    return(list(out_haz = out_haz,
                out_haz_full = out_haz_2))
}
