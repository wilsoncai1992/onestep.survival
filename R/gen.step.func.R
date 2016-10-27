#' generate step function object using y value and corresponding jump points
#'
#' @param y.vec a vector of step function values
#' @param t.vec a vector for all jump locations of step function.
#' NOTE: the first element value of y.vec is the flat part LEFT of the first jump point specified in t.vec
#'
#' @return sfun a function object, that can perform any mapping
#' @export
#'
#' @examples
#' # TO DO
gen.step.func <- function(y.vec, t.vec) {
	if (length(y.vec) != (length(t.vec) + 1)) {
		warning('the legnth of input vectors incorrect!')
	}
	sfun  <- stepfun(t.vec, y.vec, f = 0) # before the first jump point, the step function is 0 value
	return(sfun)
}
