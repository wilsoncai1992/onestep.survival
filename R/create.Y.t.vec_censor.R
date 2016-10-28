#' compute I{T.tilde >= t}
#'
#' loop over t.vec
#'
#' @param Time length n vector of failure time
#' @param t.vec t value of interest
#'
#' @return a binary vector, of length = t.vec
#' @export
#'
#' @examples
#' # TO DO
create.Y.t.vec <- function(Time, t.vec) {
	out.vec <- (Time >= t.vec) + 0
	return(out.vec)
}

#' compute I{T.tilde >= t, Delta = 1}, loop over t.vec
#'
#' @param Time length n vector of failure time
#' @param Delta length n vector of censoring indicator
#' @param t.vec t value of interest
#'
#' @return a binary vector, of length = t.vec
#' @export
#'
#' @examples
#' # TO DO
create.Y.t.vec_Delta1 <- function(Time, Delta, t.vec) {
    out.vec <- ((Time == t.vec) & (Delta == 1)) + 0
    return(out.vec)
}

