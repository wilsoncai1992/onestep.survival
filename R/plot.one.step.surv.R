#' Plot the survival curve estimator
#'
#' @param onestepfit object returned by surv.one.step or surv.one.step.complete
#' @param col line color
#' @param lty line type
#' @param add whether to add to existing plot
#'
#' @return NA
#' @export
#'
#' @examples
plot.onestep.surv <- function(onestepfit, col = 'green', lty = 1, add = FALSE) {
	if (add){
		lines(onestepfit[[1]] ~ onestepfit[[2]], col = col, lty = lty)
	}else{
		plot(onestepfit[[1]] ~ onestepfit[[2]], col = col, type = 'l', lty = lty)
	}

}
