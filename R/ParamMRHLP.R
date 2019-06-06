ParamMRHLP <- setRefClass(
  "ParamMRHLP",
  fields = list(
    fData = "FData",
    phi = "list",

    K = "numeric",
    # number of regimes
    p = "numeric",
    # dimension of beta (order of polynomial regression)
    q = "numeric",
    # dimension of w (order of logistic regression)
    variance_type = "numeric",
    nu = "numeric", # degree of freedom

    W = "matrix",
    beta = "array",
    sigma = "array"),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), K = 1, p = 2, q = 1, variance_type = 1) {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p, q = q)

      K <<- K
      p <<- p
      q <<- q
      variance_type <<- variance_type

      if (variance_type == variance_types$homoskedastic) {
        nu <<- (p + q + 3) * K - (q + 1) - (K - 1)
      }
      else{
        nu <<- (p + q + 3) * K - (q + 1)
      }

      W <<- matrix(0, q + 1, K - 1)
      beta <<- array(NA, dim = c(p + 1, fData$m, K))
      if (variance_type == variance_types$homoskedastic) {
        sigma <<- matrix(NA, fData$m, fData$m)
      }
      else{
        sigma <<- array(NA, dim = c(fData$m, fData$m, K))
      }
    },





    initParam = function(try_algo = 1) {
      "
        Initializes the regressions parameters for the MRHLP model: the
        regressions coefficients vector and the variance, for each component.
      "
      n <- nrow(phi$XBeta) # m
      m <- ncol(phi$XBeta) # P

      if (try_algo == 1) {
        #  uniform segmentation into K contiguous segments, and then a regression
        zi <- round(n / K) - 1

        s <- 0

        for (k in 1:K) {
          i <- (k - 1) * zi + 1
          j <- k * zi

          yk <- fData$Y[i:j, ]
          Xk <- phi$XBeta[i:j, ]

          beta[, , k] <<- solve(t(Xk) %*% Xk) %*% t(Xk) %*% yk

          muk <- Xk %*% beta[, , k]
          sk <- t(yk - muk) %*% (yk - muk)
          if (variance_type == variance_types$homoskedastic) {
            s <- s + sk
            sigma <<- s / n
          }
          else{
            sigma[, , k] <<- sk / length(yk)
          }
        }

      }
      else{
        # random segmentation into contiguous segments, and then a regression
        Lmin <- m + 1 #round(m/K) #nbr pts min into one segment
        tk_init <- zeros(K, 1)
        K_1 <- K
        for (k in 2:K) {
          K_1 <- K_1 - 1

          #temp <- (tk_init[k-1] + Lmin) : (m - (K_1*Lmin))

          temp <- tk_init[k - 1] + Lmin:(n - (K_1 * Lmin) - tk_init[k - 1])

          ind <- sample(length(temp))

          tk_init[k] <- temp[ind[1]]
        }
        tk_init[K + 1] <- n

        s <- 0
        for (k in 1:K) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]

          yk <- fData$Y[i:j, ]
          Xk <- phi$XBeta[i:j, ]


          beta[, , k] <<- solve(t(Xk) %*% Xk) %*% t(Xk) %*% yk

          muk <- Xk %*% beta[, , k]
          sk <- t(yk - muk) %*% (yk - muk)

          if (variance_type == variance_types$homoskedastic) {
            s <- s + sk
            sigma <<- s / n
          }
          else{
            sigma[, , k] <<- sk / length(yk)
          }
        }
      }
    },

    MStep = function(statMRHLP, verbose_IRLS) {
      # M-Step
      # Maximization w.r.t betak and sigmak (the variances)
      if (variance_type == variance_types$homoskedastic) {
        s = 0
      }
      for (k in 1:K) {
        weights <- statMRHLP$tik[, k]
        # post prob of each component k (dimension nx1)
        nk <- sum(weights)
        # expected cardinal numnber of class k

        Xk <- phi$XBeta * (sqrt(weights) %*% ones(1, p + 1))
        #[m*(p+1)]
        yk <- fData$Y * (sqrt(weights) %*% ones(1, fData$m))
        # dimension :(nx1).*(nx1) = (nx1)

        M <- t(Xk) %*% Xk
        epps <- 1e-9
        M <- M + epps * diag(p + 1)

        beta[, , k] <<- solve(M) %*% t(Xk) %*% yk # Maximization w.r.t betak
        z <- (fData$Y - phi$XBeta %*% beta[, , k]) * (sqrt(weights) %*% ones(1, fData$m))
        # Maximisation w.r.t sigmak (the variances)

        sk <- t(z) %*% z
        if (variance_type == variance_types$homoskedastic) {
          s <- s + sk

          sigma <<- s / fData$n
        } else{
          sigma[, , k] <<- sk / nk
        }
      }

      # Maximization w.r.t W
      # ----------------------------------%
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <- IRLS(phi$Xw, statMRHLP$tik, ones(nrow(statMRHLP$tik), 1), W, verbose_IRLS)

      W <<- res_irls$W
      piik <- res_irls$piik
      reg_irls <- res_irls$reg_irls
    }
  )
)
