#' @export
emMRHLP <- function(X, Y, K, p, q = 1, variance_type = 2, n_tries = 1, max_iter = 1500, threshold = 1e-6, verbose = FALSE, verbose_IRLS = FALSE) {
    fData <- FData(X, Y)

    top <- 0
    try_EM <- 0
    best_loglik <- -Inf
    cpu_time_all <- c()

    while (try_EM < n_tries) {
      try_EM <- try_EM + 1
      message("EM try nr ", try_EM)
      time <- Sys.time()

      # Initializations
      param <- ParamMRHLP$new(fData = fData, K = K, p = p, q = q, variance_type = variance_type)
      param$initParam(try_EM)
      iter <- 0
      converge <- FALSE
      prev_loglik <- -Inf

      stat <- StatMRHLP(param)

      while (!converge && (iter <= max_iter)) {
        stat$EStep(param)

        reg_irls <- param$MStep(stat, verbose_IRLS)
        stat$computeLikelihood(reg_irls)
        # FIN EM

        iter <- iter + 1
        if (verbose) {
          message("EM     : Iteration : ",
                  iter,
                  "  log-likelihood : "  ,
                  stat$log_lik)
        }
        if (prev_loglik - stat$log_lik > 1e-5) {
          message("!!!!! EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$log_lik)
          top <- top + 1
          if (top > 20)
            break
        }

        # TEST OF CONVERGENCE
        converge <-
          abs((stat$log_lik - prev_loglik) / prev_loglik) <= threshold
        if (is.na(converge)) {
          converge <-
            FALSE
        } # basicly for the first iteration when prev_loglik is Inf

        prev_loglik <- stat$log_lik
        stat$stored_loglik[iter] <- stat$log_lik
      }# FIN EM LOOP

      cpu_time_all[try_EM] <- Sys.time() - time

      # at this point we have computed param and stat that contains all the information

      if (stat$log_lik > best_loglik) {
        statSolution <- stat$copy()
        paramSolution <- param$copy()
        if (param$K == 1) {
          statSolution$tik <- matrix(stat$tik, nrow = param$fData$n, ncol = 1)
          statSolution$piik <-
            matrix(stat$piik, nrow = param$fData$n, ncol = 1)
        }
        else{
          statSolution$tik <- stat$tik[1:param$fData$n, ]
          statSolution$piik <- stat$piik[1:param$fData$n, ]
        }

        best_loglik <- stat$log_lik
      }
      if (n_tries > 1) {
        message("max value: ", stat$log_lik)
      }
    }

    # Computation of c_ig the hard partition of the curves and klas
    statSolution$MAP()

    if (n_tries > 1) {
      message("max value: ", statSolution$log_lik)
    }


    # FINISH computation of statSolution
    statSolution$computeStats(paramSolution, cpu_time_all)

    return(ModelMRHLP$new(paramMRHLP = paramSolution, statMRHLP = statSolution))
  }
