#' A Reference Class which contains parameters of a MRHLP model.
#'
#' ParamMRHLP contains all the parameters of a MRHLP model.
#'
#' @field mData [MData][MData] object representing the sample.
#' @field K The number of regimes (mixture components).
#' @field p The order of the polynomial regression.
#' @field q The dimension of the logistic regression. For the purpose of
#' segmentation, it must be set to 1.
#' @field variance_type Character indicating if the model is homoskedastic
#' (`variance_type = "homoskedastic"`) or heteroskedastic
#' (`variance_type = "heteroskedastic"`). By default the model is
#' heteroskedastic.
#' @field W Parameters of the logistic process.
#' \eqn{W = w_{1},\dots,w_{K-1}}{W = (w1,\dots,wK-1)} is a matrix of dimension
#' \eqn{(q + 1, K - 1)}, with \emph{q} the order of the logistic regression.
#' @field beta Parameters of the polynomial regressions.
#' \eqn{\beta = (\beta_{1},\dots,\beta_{K})}{\beta = (\beta1,\dots,\betaK)} is
#' a matrix of dimension \eqn{(p + 1, K)}, with \emph{p} the order of the
#' polynomial regression.
#' @field sigma2 The variances for the \emph{K} regimes. If MRHLP model is
#' homoskedastic (\emph{variance_type} = 1) then sigma2 is a matrix of size
#' \eqn{(1, 1)}, else if MRHLP model is heteroskedastic then sigma2 is a matrix
#' of size \eqn{(K, 1)}.
#' @seealso [MData]
#' @export
ParamMRHLP <- setRefClass(
  "ParamMRHLP",
  fields = list(
    mData = "MData",
    phi = "list",

    K = "numeric", # Number of regimes
    p = "numeric", # Dimension of beta (order of polynomial regression)
    q = "numeric", # Dimension of w (order of logistic regression)
    variance_type = "character",
    nu = "numeric", # Degree of freedom

    W = "matrix",
    beta = "array",
    sigma2 = "array"
  ),
  methods = list(
    initialize = function(mData = MData(numeric(1), matrix(1)), K = 1, p = 3, q = 1, variance_type = "heteroskedastic") {
      mData <<- mData

      phi <<- designmatrix(x = mData$X, p = p, q = q)

      K <<- K
      p <<- p
      q <<- q
      variance_type <<- variance_type

      if (variance_type == "homoskedastic") {
        nu <<- (q + 1) * (K - 1) + K * (p + 1) * mData$d + mData$d * (mData$d + 1) / 2
      } else {
        nu <<- (q + 1) * (K - 1) + K * (p + 1) * mData$d + K * mData$d * (mData$d + 1) / 2
      }

      W <<- matrix(0, q + 1, K - 1)
      beta <<- array(NA, dim = c(p + 1, mData$d, K))
      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA, mData$d, mData$d)
      } else {
        sigma2 <<- array(NA, dim = c(mData$d, mData$d, K))
      }
    },

    initParam = function(try_algo = 1) {
      "Method to initialize parameters \\code{W}, \\code{beta} and
      \\code{sigma2}.

      If try_algo = 1 then \\code{W}, \\code{beta} and \\code{sigma2} are
      initialized by segmenting uniformly into \\code{K} contiguous segments
      the response Y. Otherwise, \\code{W}, \\code{beta} and \\code{sigma2} are
      initialized by segmenting randomly into \\code{K} segments the response Y."

      n <- nrow(phi$XBeta)
      m <- ncol(phi$XBeta)

      if (try_algo == 1) { # Uniform segmentation into K contiguous segments, and then a regression
        zi <- round(n / K) - 1

        s <- 0

        for (k in 1:K) {
          i <- (k - 1) * zi + 1
          j <- k * zi

          yk <- mData$Y[i:j,]
          Xk <- as.matrix(phi$XBeta[i:j,])

          beta[, , k] <<- solve(t(Xk) %*% Xk) %*% t(Xk) %*% yk

          muk <- Xk %*% beta[, , k]
          sk <- t(yk - muk) %*% (yk - muk)
          if (variance_type == "homoskedastic") {
            s <- s + sk
            sigma2 <<- s / n
          } else{
            sigma2[, , k] <<- sk / length(yk)
          }
        }

      }
      else{# Random segmentation into K contiguous segments, and then a regression
        Lmin <- m + 1 # Minimum number of points in a segment
        tk_init <- zeros(K, 1)
        K_1 <- K
        for (k in 2:K) {
          K_1 <- K_1 - 1

          temp <- tk_init[k - 1] + Lmin:(n - (K_1 * Lmin) - tk_init[k - 1])

          ind <- sample(length(temp))

          tk_init[k] <- temp[ind[1]]
        }
        tk_init[K + 1] <- n

        s <- 0
        for (k in 1:K) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]

          yk <- mData$Y[i:j,]
          Xk <- phi$XBeta[i:j,]

          beta[, , k] <<- solve(t(Xk) %*% Xk) %*% t(Xk) %*% yk

          muk <- Xk %*% beta[, , k]
          sk <- t(yk - muk) %*% (yk - muk)

          if (variance_type == "homoskedastic") {
            s <- s + sk
            sigma2 <<- s / n
          }
          else{
            sigma2[, , k] <<- sk / length(yk)
          }
        }
      }
    },

    MStep = function(statMRHLP, verbose_IRLS) {
      "Method used in the EM algorithm to learn the parameters of the MRHLP model
      based on statistics provided by \\code{statMRHLP}."
      # Maximization w.r.t betak and sigmak (the variances)
      if (variance_type == "homoskedastic") {
        s = 0
      }
      for (k in 1:K) {
        weights <- statMRHLP$tik[, k] # Post probabilities of each component k (dimension nx1)
        nk <- sum(weights) # Expected cardinal numnber of class k

        Xk <- phi$XBeta * (sqrt(weights) %*% ones(1, p + 1)) # [m*(p+1)]
        yk <- mData$Y * (sqrt(weights) %*% ones(1, mData$d))

        M <- t(Xk) %*% Xk
        epps <- 1e-9
        M <- M + epps * diag(p + 1)

        beta[, , k] <<- solve(M) %*% t(Xk) %*% yk # Maximization w.r.t betak
        z <- (mData$Y - phi$XBeta %*% beta[, , k]) * (sqrt(weights) %*% ones(1, mData$d))

        # Maximisation w.r.t sigmak (the variances)
        sk <- t(z) %*% z

        if (variance_type == "homoskedastic") {
          s <- s + sk
          sigma2 <<- s / mData$m
        } else {
          sigma2[, , k] <<- sk / nk
        }
      }

      # Maximization w.r.t W
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <- IRLS(phi$Xw, statMRHLP$tik, ones(nrow(statMRHLP$tik), 1), W, verbose_IRLS)

      W <<- res_irls$W
      piik <- res_irls$piik
      reg_irls <- res_irls$reg_irls
    }
  )
)
