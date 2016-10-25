# wrapper for David's censoring SL function
# use T, A, W, Delta data format as input
# input:
# dat <- data.frame, with col: T, A, W, Delta
# T.uniq <- the unique time points of failure
# Delta.SL.Lib <- library for SuperLearner

# output:
# h.hat.t <- matrix, nrow = n.sample, ncol = length(T.uniq);
# each row is the predicted conditional hazard h(T=t | T>=t, A, W) for that subject


# 2016-09-05
# Author: Wilson
censor_SL_wrapper <- function(dat,
                              T.uniq,
                              Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth")
                              ) {
    library(plyr)
	# ----------------------------------------------------------------------------------------------------
	# transform original data into SL-friendly format
	dat_david <- dat

	# WILSON: add error catching, if there are not complete data structure
	dat_david <- rename(dat_david, c('A' = 'Z', 'T' = 'T.tilde'))

	if ('ID' %in% toupper(colnames(dat_david))) {
		# if there are already id in the dataset
		dat_david <- rename(dat_david, c('ID' = 'id'))
	}else{
		# if no id exist
		# create 'id' on our own
		dat_david$id <- 1:nrow(dat_david)
	}

	# censoring
	if ('delta' %in% colnames(dat_david)) {
	    dat_david <- rename(dat_david, c('delta' = 'Delta.J'))
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
	# perform censoring SL for the maximum time point

	# create new datalist, where censoring and outcome indicators switch each other
	datalist2 <- lapply(datalist, function(x) {x[,c('C', 'N1')] <- x[,c('N1', 'C')]; x})

	censor_hat <- estimateHazards(dataList = datalist2, J=1, verbose=FALSE, strata=NULL,
	                              adjustVars = dat_david[,baseline_name],
	                              SLlibrary.event = Delta.SL.Lib,
	                              glmFormula.event = NULL,
	                              bounds = NULL
	)
	
	# censor_hat[[1]][censor_hat[[1]]$t == T.it,'Q1Haz']
	# mean(censor_hat[[1]][censor_hat[[1]]$t == T.it,'Q1Haz'])

	# turn into wide format
	library(tidyr)
	out_censor <- censor_hat[[3]]
	out_censor <- out_censor[,c('id', 't', 'Q1Haz')]
	out_censor <- spread(out_censor, t, Q1Haz)
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