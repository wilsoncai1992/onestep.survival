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
#' @importFrom plyr rename
#' @import survtmle
#' @import tidyr
censor_SL_wrapper <- function(dat,
                              T.uniq,
                              Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")
                              ) {
	# ----------------------------------------------------------------------------------------------------
	# transform original data into SL-friendly format
	dat_david <- dat

	if ('T.TILDE' %in% toupper(colnames(dat_david))) {
	    # if there is t.tilde in variable, remove any T
	    keeps <- setdiff(colnames(dat_david), 'T')
	    dat_david <- dat_david[,keeps]
	    dat_david <- plyr::rename(dat_david, c('T.tilde' = 'ftime'))
	}else if('T' %in% toupper(colnames(dat_david))){
	    # if no t.tilde, rename to T.tilde
	    dat_david <- plyr::rename(dat_david, c('T' = 'ftime'))
	}else{
	    # if there are no T
	    stop("There should be T variable!")
	}
	dat_david <- plyr::rename(dat_david, c('A' = 'trt'))

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
	    dat_david <- plyr::rename(dat_david, c('delta' = 'ftype'))
	}else{
	    # no censoring in the dataset
	    # censoring to be all 1
	    dat_david$ftype <- 1
	}


	# remove all other useless columns
	baseline_name <- grep('W', colnames(dat_david), value = TRUE)
	keeps <- c("id", baseline_name, 'ftime', 'ftype', 'trt')
	dat_david <- dat_david[,keeps]

	# ====================================================================================================
	# remove failure time = 0
	# ====================================================================================================
	dat_david <- dat_david[dat_david$ftime != 0,]
	n.data <- nrow(dat_david)

	# ====================================================================================================
	# prepare
	# ====================================================================================================
	T.uniq <- unique(sort(dat_david$ftime))
	T.max <- max(T.uniq)

	adjustVars <- dat_david[,baseline_name]


	datalist <- survtmle::makeDataList(dat = dat_david,
	                                   J = 1, # one kind of failure
	                                   ntrt = 2, # one kind of treatment
	                                   uniqtrt = c(0,1),
	                                   t0 = T.max, # time to predict on
	                                   bounds=NULL)
	# ====================================================================================================
	# estimate g_2 (censoring)
	# perform censoring SL for the maximum time point
	# ====================================================================================================
	g2_hat <- estimateCensoring(dataList = datalist, adjustVars = adjustVars,
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
