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
plot_initial_G <- function(onestepfit, add = FALSE, ...) {
    G_hat <- colMeans(onestepfit$initial_fit$G.hat.t$out_censor_full)
    T.uniq <- onestepfit$T.uniq
    if(add){
        lines(y = G_hat, x = 1:max(T.uniq), ...)
    }else{
        plot(y = G_hat, x = 1:max(T.uniq), type = 'l', ...)
    }
}

