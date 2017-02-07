#' Title
#'
#' @param fit_obj
#' @param q
#' @param add
#' @param col
#' @param ...
#'
#' @return NA
#' @export
#'
#' @examples
plot_CI <- function(fit_obj, q = 0.95, add = FALSE, col = 'black', ...) {
    sd_CI <- sqrt(fit_obj$var)
    upper <- fit_obj$Psi.hat + qnorm(p = (1-q)/2, lower.tail = FALSE) * sd_CI
    lower <- fit_obj$Psi.hat - qnorm(p = (1-q)/2, lower.tail = FALSE) * sd_CI

    # ad-hoc thresholding between (0,1)
    upper[upper > 1] <- 1
    upper[upper < 0] <- 0
    lower[lower > 1] <- 1
    lower[lower < 0] <- 0

    step_curve_upper <- stepfun(x = fit_obj$T.uniq, y = c(1, upper))
    curve(step_curve_upper, from = 0, to = max(fit_obj$T.uniq), add = add, col = col, ...)

    step_curve_lower <- stepfun(x = fit_obj$T.uniq, y = c(1, lower))
    curve(step_curve_lower, from = 0, to = max(fit_obj$T.uniq), add = add, col = col, ...)

    return(list(upper = upper, lower = lower))
}
