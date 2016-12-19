#' Title
#'
#' @param onestepfit 
#' @param add 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
plot_initial_S <- function(onestepfit, add = FALSE, ...) {
    S_hat <- colMeans(onestepfit$initial_fit$Qn.A1.t)
    T.uniq <- onestepfit$T.uniq
    if(add){
        lines(S_hat ~ T.uniq, ...)
    }else{
        plot(S_hat ~ T.uniq, type = 'l', ...)
    }
}

