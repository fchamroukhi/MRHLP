source("R/enums.R")
source("R/utils.R")
source("R/IRLS.R")

ParamMRHLP <- setRefClass(
  "ParamMRHLP",
  fields = list(W = "matrix",
                beta = "matrix",
                sigma = "array"),
  methods = list(
    initParam = function(modelMRHLP, phi, try_algo = 1) {
      "
        Initializes the regressions parameters for the MRHLP model: the
        regressions coefficients vector and the variance, for each component.
      "
      n <- nrow(phi$XBeta) # m
      m <- ncol(phi$XBeta) # P

      if (try_algo==1){ #  uniform segmentation into K contiguous segments, and then a regression
        zi <- round(n/modelMRHLP$K)-1

        s <- 0

        for (k in 1:modelMRHLP$K){
          i <- (k-1)*zi+1
          j <- k*zi

          yk <- modelMRHLP$Y[i:j,]
          Xk <- phi$XBeta[i:j,]

          muk <- Xk %*% beta[,,k]
          sk <- t(yk - muk) %*% (yk - muk)

          if (modelMRHLP$variance_type == variance_types$homoskedastic){
            s <- s + sk
            sigma <<- s/n
          }
          else{
            sigma[, ,k] <<- sk / length(yk)
          }
        }
      }
      else{ # random segmentation into contiguous segments, and then a regression
        Lmin <- m+1 #round(m/K) #nbr pts min into one segment
        tk_init <- zeros(modelMRHLP$K,1)
        K_1 <- modelMRHLP$K
        for (k in 2:modelMRHLP$K) {
          K_1 <- K_1-1;
          #temp <- (tk_init[k-1] + Lmin) : (m - (K_1*Lmin))

          temp <- tk_init[k-1] + Lmin : (n - (K_1*Lmin) - tk_init[k-1])

          ind <- sample(length(temp));
          tk_init[k] <- temp[ind[1]]
        }
        tk_init[K+1] <- n

        s <- 0
        for (k in 1:modelMRHLP$K){
          i <- tk_init[k] + 1
          j <- tk_init[k+1]

          yk <- modelMRHLP$Y[i:j,]
          Xk <- phi$XBeta[i:j,]


          beta[,,k] <<- solve(t(Xk)%*%Xk)%*%t(Xk)%*%yk

          muk <- Xk %*% beta[,,k]
          sk <- t(yk - muk) %*% (yk - muk)

          if (modelMRHLP$variance_type == variance_types$homoskedastic){
            s <- s + sk
            sigma <<- s/n
          }
          else{
            sigma[, ,k] <<- sk / length(yk)
          }
        }
      }
    },

    MStep = function(modelMRHLP, statMRHLP, phi, verbose_IRLS) {
      # M-Step
      # Maximization w.r.t betak and sigmak (the variances)
      if (modelMRHLP$variance_type == variance_types$homoskedastic) {
        s = 0
      }
      for (k in 1:modelMRHLP$K) {
        weights <- statMRHLP$tik[, k]
        # post prob of each component k (dimension nx1)
        nk <- sum(weights)
        # expected cardinal numnber of class k

        Xk <- phi$XBeta * (sqrt(weights) %*% ones(1, modelMRHLP$p + 1))
        #[m*(p+1)]
        yk <- modelMRHLP$Y * (sqrt(weights))
        # dimension :(nx1).*(nx1) = (nx1)

        M <- t(Xk) %*% Xk
        epps <- 1e-9
        M <- M + epps * diag(modelMRHLP$p + 1)

        beta[, k] <<- solve(M) %*% t(Xk) %*% yk # Maximization w.r.t betak
        z <- sqrt(weights) * (modelMRHLP$Y - phi$XBeta %*% beta[, k])
        # Maximisation w.r.t sigmak (the variances)
        priorsigma =  0
        #1e-5;
        if (modelMRHLP$variance_type == variance_types$homoskedastic) {
          sk <- t(z) %*% z
          s <- s + sk

          sigma <<- s / modelRHLP$m
        } else{
          sigma[k] <<- t(z) %*% z / nk  + priorsigma
        }
      }

      # Maximization w.r.t W
      # ----------------------------------%
      #  IRLS : Iteratively Reweighted Least Squares (for IRLS, see the IJCNN 2009 paper)
      res_irls <-
        IRLS(statRHLP$tik,
             phi$Xw,
             W,
             verbose_IRLS = verbose_IRLS,
             piik_len = modelRHLP$m)

      W <<- res_irls$W
      piik <- res_irls$piik
      reg_irls <- res_irls$reg_irls
    }
  )
)

ParamMRHLP <- function(modelMRHLP) {
  W <- matrix(0, modelMRHLP$p + 1, modelMRHLP$K - 1)
  beta <- matrix(NA, modelMRHLP$p + 1, modelMRHLP$m, modelMRHLP$K)
  if (modelMRHLP$variance_type == variance_types$homoskedastic) {
    sigma <- matrix(NA, modelMRHLP$m, modelMRHLP$m)
  }
  else{
    sigma <- array(NA, dim = c(modelMRHLP$m, modelMRHLP$m, modelMRHLP$m))
  }
  new("ParamMRHLP", W = W, beta = beta, sigma = sigma)
}
