#' compute l2 inner product of two step functions
#'
#' f and g
#'
#' @param f.step two step-function pdf with shared jump points; can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
#' @param g.step two step-function pdf with shared jump points; can be matrix input: nrow = # of different step-function pdf, ncol = length(T.grid)
#' @param T.grid shared jump points
#'
#' @return scalar
#' @export
#'
#' @examples
#' # TO DO
l2_inner_prod_step <- function(f.step, g.step, T.grid) {
    if(is.vector(f.step) & is.vector(g.step)){
        # both f and g are one sample
        f.times.g <- f.step * g.step
    }
    if(!is.vector(f.step) & is.vector(g.step)){
        # f: multi-sample
        # g: one-sample
        f.times.g <- sweep(f.step,MARGIN=2,g.step,`*`) # multiply g to each row of f.
    }
    if(is.vector(f.step) & !is.vector(g.step)){
        # f: one-sample
        # g: multi-sample
        f.times.g <- sweep(g.step,MARGIN=2,f.step,`*`) # multiply f to each row of g.
    }
    if(!is.vector(f.step) & !is.vector(g.step)){
        # both f and g are multi-sample of same sample size
        if(nrow(f.step) != nrow(g.step)) stop('f and g have different sample size!')
        f.times.g <- f.step * g.step
    }
    # ------------------------------------------------------------------------------------
    result <- compute_step_cdf(f.times.g, T.grid)
    if(!is.vector(f.step) | !is.vector(g.step)){
        # there is multi-sample
        result <- apply(result, 1, function(obj) tail(obj, 1))
    }else{
        # both f and g are one-sample
        result <- tail(result, 1)
    }
    return(result)
}
