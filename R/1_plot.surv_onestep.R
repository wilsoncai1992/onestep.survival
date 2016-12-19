#' Plot the survival curve estimator
#'
#' @param onestepfit object returned by surv.one.step or surv.one.step.complete
#' @param col line color
#' @param add whether to add to existing plot
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
plot.surv_onestep <- function(onestepfit, col = 'green', add = FALSE, ...) {
    step_curve <- stepfun(x = onestepfit$T.uniq, y = c(1, onestepfit$Psi.hat))
    curve(step_curve, from = 0, to = max(onestepfit$T.uniq), add = add, col = col, ...)
}
