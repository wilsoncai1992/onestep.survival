#' Plot the survival curve estimator
#'
#' @param obj
#' @param add
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
plot.surv_survtmle <- function(obj, add = FALSE, ...) {
    step_curve <- stepfun(x = obj$T.uniq, y = c(1, obj$s_vec), right = TRUE)
    curve(step_curve, from = 0, to = max(obj$T.uniq), add = add, ...)
}
