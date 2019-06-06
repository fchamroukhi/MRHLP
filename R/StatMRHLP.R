#' @export
StatMRHLP <- setRefClass(
  "StatMRHLP",
  fields = list(
    piik = "matrix",
    z_ik = "matrix",
    klas = "matrix",
    Ex = "matrix",
    log_lik = "numeric",
    com_loglik = "numeric",
    stored_loglik = "list",
    stored_com_loglik = "list",
    BIC = "numeric",
    ICL = "numeric",
    AIC = "numeric",
    cpu_time = "numeric",
    log_piik_fik = "matrix",
    log_sum_piik_fik = "matrix",
    tik = "matrix",
    polynomials = "array",
    weighted_polynomials = "array"
  ),
  methods = list(
    initialize = function(paramMRHLP = ParamMRHLP(fData = FData(numeric(1), matrix(1)), K = 1, p = 2, q = 1, variance_type = 1)) {
      piik <<- matrix(NA, paramMRHLP$fData$n, paramMRHLP$K)
      z_ik <<- matrix(NA, paramMRHLP$fData$n, paramMRHLP$K)
      klas <<- matrix(NA, paramMRHLP$fData$n, 1)
      Ex <<- matrix(NA, paramMRHLP$fData$n, paramMRHLP$K)
      log_lik <<- -Inf
      com_loglik <<- -Inf
      stored_loglik <<- list()
      stored_com_loglik <<- list()
      BIC <<- -Inf
      ICL <<- -Inf
      AIC <<- -Inf
      cpu_time <<- Inf
      log_piik_fik <<- matrix(0, paramMRHLP$fData$n, paramMRHLP$K)
      log_sum_piik_fik <<- matrix(NA, paramMRHLP$fData$n, 1)
      tik <<- matrix(0, paramMRHLP$fData$n, paramMRHLP$K)
      polynomials <<- array(NA, dim = c(paramMRHLP$fData$n, paramMRHLP$p, paramMRHLP$K))
      weighted_polynomials <<- array(NA, dim = c(paramMRHLP$fData$n, paramMRHLP$p, paramMRHLP$K))
    },

    MAP = function() {
      "
      calcule une partition d'un echantillon par la regle du Maximum A Posteriori a partir des probabilites a posteriori
      Entrees : post_probas , Matrice de dimensions [n x K] des probabibiltes a posteriori (matrice de la partition floue)
      n : taille de l'echantillon
      K : nombres de classes
      klas(i) = arg   max (post_probas(i,k)) , for all i=1,...,n
      1<=k<=K
      = arg   max  p(zi=k|xi;theta)
      1<=k<=K
      = arg   max  p(zi=k;theta)p(xi|zi=k;theta)/sum{l=1}^{K}p(zi=l;theta) p(xi|zi=l;theta)
      1<=k<=K
      Sorties : classes : vecteur collones contenant les classe (1:K)
      Z : Matrice de dimension [nxK] de la partition dure : ses elements sont zik, avec zik=1 si xi
      appartient la classe k (au sens du MAP) et zero sinon.
      "
      N <- nrow(piik)
      K <- ncol(piik)
      ikmax <- max.col(piik)
      ikmax <- matrix(ikmax, ncol = 1)
      z_ik <<- ikmax %*% ones(1, K) == ones(N, 1) %*% (1:K) # partition_MAP
      klas <<- ones(N, 1)
      for (k in 1:K) {
        klas[z_ik[, k] == 1] <<- k
      }
    },
    #######
    # compute loglikelihood
    #######
    computeLikelihood = function(reg_irls) {
      log_lik <<- sum(log_sum_piik_fik) + reg_irls

    },

    #######
    # compute the final solution stats
    #######
    computeStats = function(paramMRHLP, cpu_time_all) {
      for (k in 1:K) {
        polynomials[,,k] <<- paramMRHLP$phi$XBeta %*% paramMRHLP$beta[,,k]
        weighted_polynomials[,,k] <<- (piik[,k] %*% ones(1, paramMRHLP$fData$m)) * polynomials[,,k]
      }

      #Ex <<- matrix(rowSums(weighted_polynomials))
      Ex <<- apply(weighted_polynomials, c(1,2), sum)

      cpu_time <<- mean(cpu_time_all)
      # Psi <- c(as.vector(paramRHLP$Wk), as.vector(paramRHLP$betak), as.vector(paramRHLP$sigmak))
      BIC <<- log_lik - (paramMRHLP$nu * log(paramMRHLP$fData$m) / 2)
      AIC <<- log_lik - paramMRHLP$nu


      zik_log_alphag_fg_xij <- (z_ik) * (log_piik_fik)

      com_loglik <<- sum(rowSums(zik_log_alphag_fg_xij))


      ICL <<- com_loglik - paramMRHLP$nu * log(paramMRHLP$fData$m) / 2
    },
    #######
    # EStep
    #######
    EStep = function(paramMRHLP) {
      piik <<- multinomialLogit(paramMRHLP$W, paramMRHLP$phi$Xw, ones(paramMRHLP$fData$n, paramMRHLP$K), ones(paramMRHLP$fData$n, 1))$piik
      #log_piik_fik <<- zeros(modelMRHLP$n, modelMRHLP$K)

      for (k in 1:paramMRHLP$K) {
        muk <- paramMRHLP$phi$XBeta %*% paramMRHLP$beta[,,k]
        if (variance_type == variance_types$homoskedastic) {
          sigmak <- paramMRHLP$sigma
        }else{
          sigmak <- paramMRHLP$sigma[,,k]

        }

        z <- ((paramMRHLP$fData$Y - muk) %*% solve(sigmak)) * (paramMRHLP$fData$Y - muk)

        mahalanobis <- matrix(rowSums(z))

        denom <- (2*pi)^(paramMRHLP$fData$m/2) * (det(sigmak)) ^ 0.5

        log_piik_fik[,k] <<- log(piik[,k]) - ones(paramMRHLP$fData$n, 1) %*% log(denom) - 0.5 * mahalanobis
      }

      log_piik_fik <<- pmax(log_piik_fik, log(.Machine$double.xmin))
      piik_fik <- exp(log_piik_fik)
      log_sum_piik_fik <<- matrix(log(rowSums(piik_fik)))

      log_tauik <- log_piik_fik - log_sum_piik_fik %*% ones(1,paramMRHLP$K)
      tik <<- normalize(exp(log_tauik),2)$M
    }
  )
)


