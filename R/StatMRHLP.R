#' A Reference Class which contains statistics of a MRHLP model.
#'
#' StatMRHLP contains all the parameters of a [MRHLP][ParamMRHLP] model.
#'
#' @field piik Matrix of size \eqn{(n, K)} representing the probabilities
#' \eqn{P(zi = k; W) = P(z_{ik} = 1; W)}{P(zi = k; W) = P(z_ik = 1; W)} of the
#' latent variable \eqn{zi,\ i = 1,\dots,m}{zi, i = 1,\dots,n}.
#'
#' @field z_ik Hard segmentation logical matrix of dimension \eqn{(n, K)}
#' obtained by the Maximum a posteriori (MAP) rule:
#' \eqn{z_{ik} = 1 \ \textrm{if} \ z_{ik} = \textrm{arg} \ \textrm{max}_{k} \
#' P(z_i = k | Y, W, \beta);\ 0 \ \textrm{otherwise}}{z_ik = 1 if z_ik =
#' arg max_k P(z_i = k | Y, W, \beta); 0 otherwise}, \eqn{k = 1,\dots,K}.
#' @field klas Column matrix of the labels issued from `z_ik`. Its elements are
#' \eqn{klas(i) = k}, \eqn{k = 1,\dots,K}.
#' @field Ex Column matrix of dimension \emph{n}. `Ex` is the curve expectation
#' : sum of the polynomial components \eqn{\beta_{k} \times X_{i}}{\betak x X_i}
#' weighted by the logistic probabilities `pi_ik`:
#' \eqn{Ey(i) = \sum_{k = 1}^{K} pi_ik \times \beta_{k} \times X_{i}}{Ey(i) =
#' \sum_{k=1}^K pi_ik x \betak x X_i}, \eqn{i = 1,\dots,n}.
#' @field log\_lik Numeric. Log-likelihood of the MRHLP model.
#' @field com_loglik Numeric. Complete log-likelihood of the MRHLP model.
#' @field stored_loglik List. Stored values of the log-likelihood at each EM
#' iteration.
#' @field BIC Numeric. Value of the BIC (Bayesian Information Criterion)
#' criterion. The formula is \eqn{BIC = log\_lik - nu \times
#' \textrm{log}(n) / 2}{BIC = log\_lik - nu x log(n) / 2} with \emph{nu} the
#' degree of freedom of the MRHLP model.
#' @field ICL Numeric. Value of the ICL (Integrated Completed Likelihood)
#' criterion. The formula is \eqn{ICL = com\_loglik - nu \times
#' \textrm{log}(n) / 2}{ICL = com_loglik - nu x log(n) / 2} with \emph{nu} the
#' degree of freedom of the MRHLP model.
#' @field AIC Numeric. Value of the AIC (Akaike Information Criterion)
#' criterion. The formula is \eqn{AIC = log\_lik - nu}{AIC = log\_lik - nu}.
#' @field cpu_time Numeric. Average executon time of a EM step.
#' @field log_piik_fik Matrix of size \eqn{(n, K)} giving the values of the
#' logarithm of the joint probability
#' \eqn{P(Y_{i}, \ zi = k)}{P(Yi, zi = k)}, \eqn{i = 1,\dots,n}.
#' @field log_sum_piik_fik Column matrix of size \emph{n} giving the values of
#' \eqn{\sum_{k = 1}^{K} \textrm{log} P(Y_{i}, \ zi = k)}{\sum_{k = 1}^{K} log
#' P(Yi, zi = k)}, \eqn{i = 1,\dots,n}.
#' @field tik Matrix of size \eqn{(n, K)} giving the posterior probability that
#' \eqn{Y_{i}}{Yi} originates from the \eqn{k}-th regression model
#' \eqn{P(zi = k | Y, W, \beta)}.
#' @field polynomials Matrix of size \eqn{(n, K)} giving the values of
#' \eqn{\beta_{k} \times X_{i}}{\betak x X_i}, \eqn{i = 1,\dots,n}.
#' @seealso [ParamMRHLP], [FData]
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

    initialize = function(paramMRHLP = ParamMRHLP()) {

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
      "Method used in the EM algorithm to update statistics based on parameters
      provided by \\code{paramMRHLP} (prior and posterior probabilities)."
      piik <<- multinomialLogit(paramMRHLP$W, paramMRHLP$phi$Xw, ones(paramMRHLP$fData$n, paramMRHLP$K), ones(paramMRHLP$fData$n, 1))$piik
      #log_piik_fik <<- zeros(modelMRHLP$n, modelMRHLP$K)

      for (k in 1:paramMRHLP$K) {
        muk <- paramMRHLP$phi$XBeta %*% paramMRHLP$beta[,,k]
        if (variance_type == "homoskedastic") {
          sigma2k <- paramMRHLP$sigma2
        }else{
          sigma2k <- paramMRHLP$sigma2[,,k]

        }

        z <- ((paramMRHLP$fData$Y - muk) %*% solve(sigma2k)) * (paramMRHLP$fData$Y - muk)

        mahalanobis <- matrix(rowSums(z))

        denom <- (2*pi)^(paramMRHLP$fData$m/2) * (det(sigma2k)) ^ 0.5

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


