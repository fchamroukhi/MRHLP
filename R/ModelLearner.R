#' emMRHLP is used to fit a MRHLP model.
#'
#' emMRHLP is used to fit a MRHLP model. The estimation method is performed by
#' the Expectation-Maximization algorithm.
#'
#' @details emMRHLP function is based on the EM algorithm. This functions starts
#' with an initialization of the parameters done by the method `initParam` of
#' the class [ParamMRHLP][ParamMRHLP], then it alternates between a E-Step
#' (method of the class [StatMRHLP][StatMRHLP]) and a M-Step (method of the class
#' [ParamMRHLP][ParamMRHLP]) until convergence (until the absolute difference of
#' log-likelihood between two steps of the EM algorithm is less than the
#' `threshold` parameter).
#'
#' @param X Numeric vector of length \emph{m} representing the covariates.
#' @param Y Matrix of size \eqn{(n, m)} representing \emph{n} functions of `X`
#' observed at points \eqn{1,\dots,m}.
#' @param K The number of regimes (mixture components).
#' @param p The order of the polynomial regression.
#' @param q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @param variance_type Optional character indicating if the model is
#' "homoskedastic" or "heteroskedastic". By default the model is
#' "heteroskedastic".
#' @param n_tries Number of times EM algorithm will be launched.
#' The solution providing the highest log-likelihood will be returned.
#'
#' If `n_tries` > 1, then for the first pass, parameters are initialized
#' by uniformly segmenting the data into K segments, and for the next passes,
#' parameters are initialized by randomly segmenting the data into K contiguous
#'  segments.
#' @param max_iter The maximum number of iterations for the EM algorithm.
#' @param threshold A numeric value specifying the threshold for the relative
#'  difference of log-likelihood between two steps  of the EM as stopping
#'  criteria.
#' @param verbose A logical value indicating whether values of the
#' log-likelihood should be printed during EM iterations.
#' @param verbose_IRLS A logical value indicating whether values of the
#' criterion optimized by IRLS should be printed at each step of the EM
#' algorithm.
#' @return EM returns an object of class [ModelMRHLP][ModelMRHLP].
#' @seealso [ModelMRHLP], [ParamMRHLP], [StatMRHLP]
#' @export
emMRHLP <- function(X, Y, K, p = 3, q = 1, variance_type = c("heteroskedastic", "homoskedastic"), n_tries = 1, max_iter = 1500, threshold = 1e-6, verbose = FALSE, verbose_IRLS = FALSE) {
    if (is.vector(Y)) { # Univariate time series
      Y <- as.matrix(Y)
    }
    mData <- MData(X, Y)

    top <- 0
    try_EM <- 0
    best_loglik <- -Inf
    cpu_time_all <- c()

    while (try_EM < n_tries) {
      try_EM <- try_EM + 1

      if (verbose) {
        cat(paste0("EM try number: ", try_EM, "\n\n"))
      }

      time <- Sys.time()

      # Initialization
      variance_type <- match.arg(variance_type)
      param <- ParamMRHLP$new(mData = mData, K = K, p = p, q = q, variance_type = variance_type)
      param$initParam(try_EM)
      iter <- 0
      converge <- FALSE
      prev_loglik <- -Inf

      stat <- StatMRHLP(param)

      while (!converge && (iter <= max_iter)) {
        stat$EStep(param)

        reg_irls <- param$MStep(stat, verbose_IRLS)
        stat$computeLikelihood(reg_irls)

        iter <- iter + 1
        if (verbose) {
          cat(paste0("EM: Iteration : ", iter, " || log-likelihood : ", stat$log_lik, "\n"))
        }
        if (prev_loglik - stat$log_lik > 1e-5) {
          warning(paste0("EM log-likelihood is decreasing from ", prev_loglik, "to ", stat$log_lik, " !"))
          top <- top + 1
          if (top > 20)
            break
        }

        # Test of convergence
        converge <- abs((stat$log_lik - prev_loglik) / prev_loglik) <= threshold
        if (is.na(converge)) {
          converge <-
            FALSE
        } # Basically for the first iteration when prev_loglik is Inf

        prev_loglik <- stat$log_lik
        stat$stored_loglik[iter] <- stat$log_lik
      } # End of the EM loop

      cpu_time_all[try_EM] <- Sys.time() - time

      if (stat$log_lik > best_loglik) {
        statSolution <- stat$copy()
        paramSolution <- param$copy()
        if (param$K == 1) {
          statSolution$tik <- matrix(stat$tik, nrow = param$mData$m, ncol = 1)
          statSolution$piik <- matrix(stat$piik, nrow = param$mData$m, ncol = 1)
        } else {
          statSolution$tik <- stat$tik[1:param$mData$m,]
          statSolution$piik <- stat$piik[1:param$mData$m,]
        }

        best_loglik <- stat$log_lik
      }
      if (n_tries > 1) {
        cat(paste0("Max value of the log-likelihood: ", stat$log_lik, "\n"))
      }
    }

    # Computation of c_ig the hard partition of the curves and klas
    statSolution$MAP()

    if (n_tries > 1) {
      cat(paste0("Max value of the log-likelihood: ", statSolution$log_lik, "\n"))
    }


    # Finish the computation of statistics
    statSolution$computeStats(paramSolution, cpu_time_all)

    return(ModelMRHLP$new(paramMRHLP = paramSolution, statMRHLP = statSolution))
  }
