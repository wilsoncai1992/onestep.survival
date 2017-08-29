#' One-step TMLE estimator for survival curve
#'
#' one-step TMLE estimate of the treatment specific survival curve. Under right-censored data
#'
#' options to ADD:
#' SL.formula: the covariates to include in SL
#'
#' @param dat A data.frame with columns T.tilde, delta, A, W. T.tilde = min(T, C) is either the failure time of censor time, whichever happens first. 'delta'= I(T <= C) is the indicator of whether we observe failure time. A is binary treatment. W is baseline covariates. All columns with character "W" will be treated as baseline covariates.
#' @param dW A binary vector specifying dynamic treatment (as a function output of W)
#' @param g.SL.Lib A vector of string. SuperLearner library for fitting treatment regression
#' @param Delta.SL.Lib A vector of string. SuperLearner library for fitting censoring regression
#' @param ht.SL.Lib A vector of string. SuperLearner library for fitting conditional hazard regression
#' @param epsilon.step numeric. step size for one-step recursion
#' @param max.iter integer. maximal number of recursion for one-step
#' @param tol numeric. tolerance for optimization
#' @param T.cutoff int. Enforce randomized right-censoring to the observed data, so that don't estimate survival curve beyond a time point. Useful when time horizon is long.
#' @param verbose boolean. When TRUE, plot the initial fit curve, and output the objective function value during optimzation
#' @param ... additional options for plotting initial fit curve
#'
#' @return Psi.hat A numeric vector of estimated treatment-specific survival curve
#' @return T.uniq A vector of descrete time points where Psi.hat take values (have same length as Psi.hat)
#' @return params A list of estimation parameters set by user
#' @return variables A list of data summary statistics
#' @return initial_fit A list of initial fit values (hazard, g_1, Delta)
#'
#' @export
#'
#' @examples
#' library(simcausal)
#' D <- DAG.empty()
#'
#' D <- D +
#'     node("W", distr = "rbinom", size = 1, prob = .5) +
#'     node("A", distr = "rbinom", size = 1, prob = .15 + .5*W) +
#'     node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
#'     node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
#'     node("T", distr = "rconst", const = round(Trexp*100,0)) +
#'     node("C", distr = "rconst", const = round(Cweib*100, 0)) +
#'     node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
#'     node("delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
#' setD <- set.DAG(D)
#'
#' dat <- sim(setD, n=3e2)
#'
#' library(dplyr)
#' # only grab ID, W's, A, T.tilde, Delta
#' Wname <- grep('W', colnames(dat), value = TRUE)
#' dat <- dat[,c('ID', Wname, 'A', "T.tilde", "delta")]
#'
#' dW <- rep(1, nrow(dat))
#' onestepfit <- surv_onestep(dat = dat,
#'                             dW = dW,
#'                             verbose = FALSE,
#'                             epsilon.step = 1e-3,
#'                             max.iter = 1e3)
#' @import dplyr
#' @import survtmle
#' @import abind
#' @import SuperLearner
surv_onestep <- function(dat,
                          dW = rep(1, nrow(dat)),
                          g.SL.Lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                          Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          epsilon.step = 1e-5,
                          max.iter = 1e3,
                          tol = 1/nrow(dat),
                          T.cutoff = NULL,
                          verbose = TRUE,
                          ...) {
    # ===================================================================================
    # preparation
    # ===================================================================================
    after_check <- check_and_preprocess_data(dat = dat, dW = dW, T.cutoff = T.cutoff)
    dat <- after_check$dat
    dW <- after_check$dW
    n.data <- after_check$n.data
    W_names <- after_check$W_names

    W <- dat[,W_names]
    W <- as.data.frame(W)

    # dW check
    if(all(dW == 0)) {
        dat$A <- 1 - dat$A # when dW is all zero
        dW <- 1 - dW
    }else if(all(dW == 1)){

    }else{
        stop('not implemented!')
    }

    T.uniq <- sort(unique(dat$T.tilde))
    T.max <- max(T.uniq)
    # ===================================================================================
    # estimate g(A|W)
    # ===================================================================================
    gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=g.SL.Lib, family="binomial")
    # g.hat for each observation
    g.fitted <- gHatSL$SL.predict
    # ===================================================================================
    # conditional hazard (by SL)
    # ===================================================================================
    message('estimating conditional hazard')

    h.hat.t <- estimate_hazard_SL(dat = dat, T.uniq = T.uniq, ht.SL.Lib = ht.SL.Lib)
    # h.hat at all time t=[0,t.max]
    h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
    # h.hat at observed unique time t = T.grid
    h.hat.t <- as.matrix(h.hat.t$out_haz)
    # ===================================================================================
    # estimate censoring G(A|W)
    # ===================================================================================
    message('estimating censoring')
    G.hat.t <- estimate_censoring_SL(dat = dat, T.uniq = T.uniq,
                                 Delta.SL.Lib = Delta.SL.Lib)
    # cutoff <- 0.1
    cutoff <- 0.05
    if(any(G.hat.t$out_censor_full <= cutoff)){
        warning('G.hat has extreme small values! lower truncate to 0.05')
        G.hat.t$out_censor_full[G.hat.t$out_censor_full < cutoff] <- cutoff
        G.hat.t$out_censor[G.hat.t$out_censor < cutoff] <- cutoff
    }

    Gn.A1.t_full <- as.matrix(G.hat.t$out_censor_full)
    Gn.A1.t <- as.matrix(G.hat.t$out_censor)
    # ===================================================================================
    # Gn.A1.t
    # ===================================================================================
    # plot initial fit
    if (verbose) lines(colMeans(Gn.A1.t) ~ T.uniq, col = 'yellow', lty = 1)
    # ===================================================================================
    # Qn.A1.t
    # ===================================================================================
    Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

    # ------------------------------------------------------------------------------------
    # compute cumulative hazard
    # ------------------------------------------------------------------------------------
    # cum-product approach (2016-10-05)
    Qn.A1.t_full <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full))
    for (it in 1:n.data) {
        Qn.A1.t_full[it,] <- cumprod(1 - h.hat.t_full[it,])
    }
    Qn.A1.t <- Qn.A1.t_full[,T.uniq]

    # plot initial fit
    if (verbose) lines(colMeans(Qn.A1.t) ~ T.uniq, ...)
    # ===================================================================================
    # qn.A1.t
    # ===================================================================================
    # WILSON: rewrite in sweep?
    qn.A1.t_full <- matrix(0, nrow = n.data, ncol = ncol(Qn.A1.t_full))
    for (it.n in 1:n.data) {
        qn.A1.t_full[it.n,] <- h.hat.t_full[it.n,] * Qn.A1.t_full[it.n,]
    }
    qn.A1.t <- qn.A1.t_full[,T.uniq]

    # ===================================================================================
    # D1.t: calculate IC
    # D1.A1.t: calculate IC under intervention
    # ===================================================================================

    compute_IC <- function(dat, dW, T.uniq, h.hat.t_full, g.fitted, Gn.A1.t_full, Qn.A1.t, Qn.A1.t_full) {
        I.A.dW <- dat$A == dW
        n.data <- nrow(dat)

        D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
        D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

        for (it.n in 1:n.data) {

            t_Delta1.vec <- create_Yt_vector_with_censor(Time = dat$T.tilde[it.n], Delta = dat$delta[it.n], t.vec = 1:max(T.uniq))
            t.vec <- create_Yt_vector(Time = dat$T.tilde[it.n], t.vec = 1:max(T.uniq))
            alpha2 <- (t_Delta1.vec - t.vec * h.hat.t_full[it.n,])

            alpha1 <- -I.A.dW[it.n]/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]
            alpha1_A1 <- -1/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]

            not_complete <- alpha1 * alpha2
            not_complete_A1 <- alpha1_A1 * alpha2
            # D1 matrix
            D1.t[it.n, ] <- cumsum(not_complete)[T.uniq] * Qn.A1.t[it.n,] # complete influence curve
            D1.A1.t[it.n, ] <- cumsum(not_complete_A1)[T.uniq] * Qn.A1.t[it.n,] # also update those A = 0.
        }

        # turn unstable results to 0
        D1.t[is.na(D1.t)] <- 0
        D1.A1.t[is.na(D1.A1.t)] <- 0

        return(list(D1.t = D1.t,
                    D1.A1.t = D1.A1.t))
    }

    initial_IC <- compute_IC(dat = dat,
                             dW = dW,
                             T.uniq = T.uniq,
                             h.hat.t_full = h.hat.t_full,
                             g.fitted = g.fitted,
                             Gn.A1.t_full = Gn.A1.t_full,
                             Qn.A1.t = Qn.A1.t,
                             Qn.A1.t_full = Qn.A1.t_full)
    D1.t <- initial_IC$D1.t
    D1.A1.t <- initial_IC$D1.A1.t
    # ===================================================================================
    # Pn.D1: efficient IC average
    # ===================================================================================
    # Pn.D1 vector
    Pn.D1.t <- colMeans(D1.t)
    # ===================================================================================
    # update
    # ===================================================================================
    message('targeting')
    stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17
    if(verbose) print(stopping.criteria)

    update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
    iter.count <- 0
    stopping.prev <- Inf
    all_stopping <- numeric(stopping.criteria)
    all_loglikeli <- numeric()

    while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
    # while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
        if(verbose) print(stopping.criteria)
        # =============================================================================
        # update the qn
        # ------------------------------------------------------------------------
        # vectorized
        update.mat <- compute_onestep_update_matrix(D1.t.func.prev = D1.t,
                                     Pn.D1.func.prev = Pn.D1.t,
                                     dat = dat,
                                     T.uniq = T.uniq,
                                     W_names = W_names,
                                     dW = dW)
        # ------------------------------------------------------------------------
        update.tensor <- update.tensor + update.mat

        # accelerate when log-like becomes flat
        # if((stopping.prev - stopping.criteria) > 0 & (stopping.prev - stopping.criteria) < 1e-3) update.tensor <- update.tensor + update.mat*10

        # intergrand <- rowSums(update.tensor)
        # intergrand <- apply(update.tensor, c(1,2), sum)
        intergrand <- update.tensor
        intergrand[is.na(intergrand)] <- 0
        qn.current <- qn.A1.t * exp(epsilon.step * intergrand)
        qn.current_full <- qn.A1.t_full * exp(epsilon.step * replicate(T.max, intergrand[,1])) #10-23

        # For density sum > 1: normalize the updated qn
        norm.factor <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf)[,1] #09-06
        qn.current[norm.factor > 1,] <- qn.current[norm.factor > 1,] / norm.factor[norm.factor > 1] #09-06
        qn.current_full[norm.factor > 1,] <- qn.current_full[norm.factor > 1,] / norm.factor[norm.factor > 1] #10-23

        # 11-26
        # For density sum > 1: truncate the density outside sum = 1 to be zero
        # i.e. flat cdf beyond sum to 1
        # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = -Inf)
        # qn.current[cdf_per_subj > 1] <- 0
        # cdf_per_subj <- compute_step_cdf(pdf.mat = qn.current_full, t.vec = 1:max(T.uniq), start = -Inf)
        # qn.current_full[cdf_per_subj > 1] <- 0

        # if some qn becomes all zero, prevent NA exisitence
        qn.current[is.na(qn.current)] <- 0
        qn.current_full[is.na(qn.current_full)] <- 0 #10-23
        # =============================================================================
        # compute new Qn

        Qn.current <- compute_step_cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf) # 2016-09-06
        cdf_offset <- 1 - Qn.current[,1] # 2016-09-06
        Qn.current <- Qn.current + cdf_offset # 2016-09-06

        Qn.current_full <- compute_step_cdf(pdf.mat = qn.current_full, t.vec = 1:max(T.uniq), start = Inf) # 10-23
        cdf_offset <- 1 - Qn.current_full[,1] # 10-23
        Qn.current_full <- Qn.current_full + cdf_offset # 10-23

        # check error
        # all.equal(compute_step_cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
        # ------------------------------------------------
        # print('New Qn')
        # =============================================================================
        # compute new h_t
        h.hat.t_full_current <- matrix(0, nrow = n.data, ncol = max(T.uniq))
        for (it.n in 1:n.data) {
            h.hat.t_full_current[it.n, ] <- qn.current_full[it.n, ] / Qn.current_full[it.n,]
        }
        # ------------------------------------------------
        # print('New D1')
        # compute new D1
        updated_IC <- compute_IC(dat = dat,
                                 dW = dW,
                                 T.uniq = T.uniq,
                                 h.hat.t_full = h.hat.t_full_current,
                                 g.fitted = g.fitted,
                                 Gn.A1.t_full = Gn.A1.t_full,
                                 Qn.A1.t = Qn.current,
                                 Qn.A1.t_full = Qn.current_full)

        D1.t <- updated_IC$D1.t
        D1.A1.t <- updated_IC$D1.A1.t

        # compute new Pn.D1
        Pn.D1.t <- colMeans(D1.t)
        # ===================================================================================
        # previous stopping criteria
        stopping.prev <- stopping.criteria
        # new stopping criteria
        # stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq)
        stopping.criteria <- sqrt(l2_inner_prod_step(Pn.D1.t, Pn.D1.t, T.uniq)/max(T.uniq))
        iter.count <- iter.count + 1
        # ===================================================================================
        # evaluate log-likelihood

        # construct obj
        obj <- list()
        obj$qn.current_full <- qn.current_full
        obj$Qn.current_full <- Qn.current_full
        obj$h.hat.t_full_current <- h.hat.t_full_current
        obj$dat <- dat
        # eval loglikeli
        loglike_here <- eval_loglike(obj, dW)
        all_loglikeli <- c(all_loglikeli, loglike_here)
        all_stopping <- c(all_stopping, stopping.criteria)


        ########################################################################
        # FOR DEBUG ONLY
        # if (TRUE) {
        #   # ------------------------------------------------------------
        #   # q <- seq(0,10,.1)
        #   ## truesurvExp <- 1 - pexp(q, rate = 1)
        #   # truesurvExp <- 1 - pexp(q, rate = .5)
        #   # plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red', main = paste('l2 error =', stopping.criteria))
        #
        #   # library(survival)
        #   # n.data <- nrow(dat)
        #   # km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
        #   # lines(km.fit)
        #   # ------------------------------------------------------------
        #   # Psi.hat <- colMeans(Qn.current[dat$A==dW,])
        #
        #   # Psi.hat <- colMeans(Qn.current[dat$A==dW & dat$W==0,])
        #   # ------------------------------------------------------------
        #   # Q_weighted <- Qn.current/g.fitted[,1] # 09-18: inverse weight by propensity score
        #   # Q_weighted[dat$A!=dW,] <- 0
        #   # Psi.hat <- colMeans(Q_weighted) # 09-18: inverse weight by propensity score
        #   # ------------------------------------------------------------
        #     # 10-06: update all subjects with same W strata
        #   Psi.hat <- colMeans(Qn.current)
        #   # ------------------------------------------------------------
        #
        #   lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
        #   # ------------------------------------------------------------
        #   # legend('topright', lty=1, legend = c('true', 'KM', 'one-step'), col=c('red', 'black', 'blue'))
        # }
        ########################################################################
        if (iter.count == max.iter) {
            warning('Max Iter count reached, stop iteration.')
        }
    }

    if (!exists('Qn.current')) {
        # if the iteration immediately converge
        message('converge suddenly!')
        Qn.current <- Qn.A1.t
        updated_IC <- initial_IC
        Psi.hat <- colMeans(Qn.current)
    }

    # ===================================================================================
    # compute the target parameter
    # ===================================================================================
    # return the mean of those with observed A == dW
    Psi.hat <- colMeans(Qn.current)
    # variance of the EIC
    var_CI <- apply(updated_IC$D1.t, 2, var)/n.data
    # --------------------------------------------------
    # sup norm for each dim of EIC
    sup_norm_EIC <- abs(Pn.D1.t)

    variables <- list(T.uniq = T.uniq, Qn.current = Qn.current, D1.A1.t = D1.A1.t, D1.t = D1.t, Pn.D1.t = Pn.D1.t, sup_norm_EIC = sup_norm_EIC)
    params <- list(stopping.criteria = stopping.criteria, epsilon.step = epsilon.step, iter.count = iter.count, max.iter = max.iter, dat = dat, dW = dW)
    initial_fit <- list(h.hat.t = h.hat.t, Qn.A1.t = Qn.A1.t, qn.A1.t = qn.A1.t, G.hat.t = G.hat.t, g.fitted = g.fitted)
    convergence <- list(all_loglikeli = all_loglikeli, all_stopping = all_stopping)
    # --------------------------------------------------
    to.return <- list(Psi.hat = Psi.hat,
                      T.uniq = T.uniq,
                      var = var_CI,
                      params = params,
                      variables = variables,
                      initial_fit = initial_fit,
                      convergence = convergence)
    class(to.return) <- 'surv_onestep'
    return(to.return)

}
