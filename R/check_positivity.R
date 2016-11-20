#' Check positivity assumption of input data
#'
#' Check whether positivity assumption is violated in data
#' censored data NOT supported
#'
#' @param dat input data.frame
#' @param posit_level level of positivity to detect from data
#'
#' @return raise warning when positivity violated
#' @export
#'
#' @examples
check_positivity <- function(dat, posit_level = 0.05) {
	W_names <- grep('W', colnames(dat), value = TRUE)
	W <- dat[,W_names]
	A <- dat[,'A']

	table_out <- ftable(cbind(W,A))
	if (is.null(dim(W))) {
		# if W is vector
		table_out <- table(data.frame(cbind( W,A)))
	}
	table_out <- table_out/rowSums(table_out)
	# there are strata where there is NO observation
	table_out <- table_out[complete.cases(table_out),]
	# table_out[is.na(table_out)] <- 0

	posit_violate_by_strata <- rowSums(table_out < posit_level)
	# some strata have absolutely no entry
	# posit_violate_by_strata[posit_violate_by_strata == 2] <- 1

	if(any(as.logical(posit_violate_by_strata))){
		warning(paste('Positivity assumption violated', sum(posit_violate_by_strata), 'out of', length(posit_violate_by_strata), 'stratas!'))
	}
}
