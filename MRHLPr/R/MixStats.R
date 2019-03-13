MixStats <- setRefClass(
  "MixStats",
  fields = list(
    h_ig="matrix", # h_ig = prob(curve|cluster_g) : post prob (fuzzy segmentation matrix of dim [nxG])
    c_ig="matrix", # c_ig : Hard partition obtained by the AP rule :  c_{ig} = 1
                   # if and only c_i = arg max_g h_ig (g=1,...,G)
    klas="matrix", # klas : column vector of cluster labels
    Ex_g = "matrix", # Ex_g: curve expectation: sum of the polynomial components beta_gk ri weighted by
                     # the logitic probabilities pij_gk: Ex_g(j) = sum_{k=1}^K pi_jgk beta_gk rj, j=1,...,m. Ex_g
                     # is a column vector of dimension m for each g.
    log_lik="numeric",  # the loglikelihood of the EM or CEM algorithm
    com_loglik="numeric", # the complete loglikelihood of the EM (computed at the convergence) or CEM algorithm
    stored_loglik = "list", # vector of stored valued of the comp-log-lik at each EM teration
    stored_com_loglik = "list",
    tau_ijgk = "array", # tau_ijgk prob(y_{ij}|kth_segment,cluster_g), fuzzy
                        # segmentation for the cluster g. matrix of dimension
                        # [nmxK] for each g  (g=1,...,G).
    log_tau_ijgk = "array",
    BIC="numeric", # BIC value = loglik - nu*log(nm)/2.
    ICL="numeric", # ICL value = comp-loglik_star - nu*log(nm)/2.
    AIC="numeric", # AIC value = loglik - nu.
    cpu_time = "numeric",
    log_fg_xij="matrix",
    log_alphag_fg_xij="matrix",
    polynomials="array",
    weighted_polynomials="array"
  ),
  methods=list(
    showDataClusterSegmentation = function(mixModel, mixParamSolution){

      t <- seq(0,mixModel$m-1)
      colors = c('red','blue','green', 'aquamarine', 'bisque3', 'cyan', 'darkorange', 'darkorchid', 'gold')
      colors_cluster_means = c('red','blue','green', 'aquamarine', 'bisque3', 'cyan', 'darkorange', 'darkorchid', 'gold')

      par(mfrow=c(1,1))
      for (g in 1:mixModel$G){
        cluster_g = mixModel$X[klas==g,]

        if (g==1){
          plot(t, t(cluster_g)[,1], type = 'l', col=colors[g], xlab='Time', ylab = 'y')
        }
        else{
          lines(t, t(cluster_g)[,1], type = 'l', col=colors[g])
        }

        for (i in 1:nrow(cluster_g)){
          #print(colors[i])
          lines(t, t(cluster_g)[,i], col=colors[g])
        }

      }

      for(g in 1:mixModel$G){
        lines(t,Ex_g[,g], col=colors_cluster_means[g], lwd=5)
      }


      ##########################################################################"
      ##########################################################################
      ##########################################################################

      for (g in 1:mixModel$G){
        par(mfrow=c(2,1))
        cluster_g = mixModel$X[klas==g,]
        plot(t, t(cluster_g)[,1], type = 'l', col=colors[g], ylab = 'y')
        for (k in 1:K){
          lines(t,polynomials[,k,g], lty =2 , col="black",lwd=1)
        }
        lines(t,Ex_g[,g], col=colors_cluster_means[g], lwd=5)

        for (k in 1:mixModel$K){
          if (k==1){
            plot(t,mixParamSolution$pi_jgk[1:mixModel$m,k,g], xlab='Time', ylab = 'Logistic proportions', type = "l", lwd=2, col=colors_cluster_means[k])
          }
          else{
            lines(t,mixParamSolution$pi_jgk[1:mixModel$m,k,g] , lwd=2, col=colors_cluster_means[k])
          }
        }
      }


    },

    MAP = function(){
      "
         calculate a partition by applying the Maximum A Posteriori Bayes
         allocation rule
      "
      N <- nrow(h_ig)
      K <- ncol(h_ig)
      ikmax <- max.col(h_ig)
      ikmax <- matrix(ikmax, ncol = 1)
      c_ig <<- ikmax%*%ones(1,K) == ones(N,1)%*%(1:K)
      klas <<- ones(N,1)
      for (k in 1:K){
        klas[c_ig[,k]==1] <<- k
      }
    },

    #######
    # compute the final solution stats
    #######
    computeStats = function(mixModel, mixParam, phi, cpu_time_all){
      for (g in 1:mixModel$G){
        polynomials[,,g] <<- phi$XBeta[1:mixModel$m, ] %*% mixParam$betag[,,g]
        if (K!=1 && G!=1){
          weighted_polynomials[,,g] <<- mixParam$pi_jgk[,,g] * polynomials[,,g]
          Ex_g[,g] <<- rowSums(weighted_polynomials[,,g])
        }
        else if (K==1 && G!=1){
          weighted_polynomials[,,g] <<- mixParam$pi_jgk[,g] * polynomials[,,g]
          Ex_g[,g] <<- weighted_polynomials[,,g]
        }
        else if (K!=1 && G==1){
          weighted_polynomials[,,g] <<- mixParam$pi_jgk * polynomials[,,g]
          Ex_g[,g] <<- matrix(rowSums(weighted_polynomials[,,g]))
        }
        else{ #(K==1 && G==1)
          weighted_polynomials[,,g] <<- mixParam$pi_jgk * polynomials[,,g]
          Ex_g[,g] <<- weighted_polynomials[,,g]
        }
      }

      Ex_g <<- matrix(Ex_g, nrow = mixModel$m)
      cpu_time <<- mean(cpu_time_all)
      Psi <- c(as.vector(mixParam$alpha_g), as.vector(mixParam$Wg), as.vector(mixParam$betag), as.vector(mixParam$sigmag))
      nu <- length(Psi)
      BIC <<- log_lik - (nu*log(mixModel$n)/2)
      AIC <<- log_lik - nu

      cig_log_alphag_fg_xij <- (c_ig)*(log_alphag_fg_xij);
      com_loglik <<- sum(rowSums(cig_log_alphag_fg_xij));

      ICL <<- com_loglik - nu*log(mixModel$n)/2;

    },

    #######
    # EStep
    #######

    EStep = function(mixModel, mixParam, phi, variance_type){
      piik <- modele_logit(mixParam$Wg, phi$Xw)

      log_piik_fik <- zeros(mixModel$m, mixModel$K)

      for (k in 1:mixModel$K){
        muk <- phi$XBeta %*% mixParam$betag[,,k]
        if (variance_type == variance_types.homoskedastic){
          sigmak <- mixParam$sigmag
        }else{
          sigmak <- mixParam$sigmag[,,k]
        }

        z <- ((mixModel$y - muk) %*% solve(sigmak)) * (mixModel$y - muk)

        mahalanobis <- rowSums(z)

        denom <- (2*pi)^(mixModel$d/2) * (det(sigmak)) ^ 0.5

        log_piik_fik[,k] <- log(piik[,k]) - ones(mixModel$m, 1) %*% log(denom) - 0.5* mahalanobis
      }

      log_piik_fik <<- pmax(log_piik_fik, log(.Machine$double.xmin))
      piik_fik <- exp(log_piik_fik)
      log_sum_piik_fik <- log(rowSum(piik_fik))
      log_tauik <- log_piik_fik - log_sum_piik_fik %*% ones(1,mixModel$K)
      tauik <<- normalize(exp(log_tauik),2)

    }

  )
)


MixStats<-function(mixModel, options){
  h_ig <- matrix(NA,mixModel$n, mixModel$G)
  c_ig <- matrix(NA,mixModel$n, mixModel$G)
  klas <- matrix(NA, mixModel$n,1)
  Ex_g <- matrix(NA,mixModel$m, mixModel$G)
  log_lik <- -Inf
  com_loglik <- -Inf
  stored_loglik <- list()
  stored_com_loglik <- list()
  BIC <- -Inf
  ICL <- -Inf
  AIC <- -Inf
  cpu_time <- Inf
  log_fg_xij <- matrix(0, mixModel$n, mixModel$G)
  log_alphag_fg_xij <- matrix(0, mixModel$n, mixModel$G)
  polynomials <- array(NA, dim = c(mixModel$m, mixModel$K, mixModel$G))
  weighted_polynomials <- array(NA, dim = c(mixModel$m, mixModel$K, mixModel$G))
  tau_ijgk <- array(0, dim = c(mixModel$n*mixModel$m, mixModel$K, mixModel$G))
  log_tau_ijgk <- array(0, dim = c(mixModel$n*mixModel$m, mixModel$K, mixModel$G))

  new("MixStats", h_ig=h_ig, c_ig=c_ig, klas=klas, Ex_g=Ex_g, log_lik=log_lik, com_loglik=com_loglik, stored_loglik=stored_loglik, stored_com_loglik=stored_com_loglik, BIC=BIC, ICL=ICL, AIC=AIC, cpu_time=cpu_time,
      log_fg_xij=log_fg_xij, log_alphag_fg_xij=log_alphag_fg_xij, polynomials=polynomials, weighted_polynomials=weighted_polynomials, tau_ijgk=tau_ijgk, log_tau_ijgk=log_tau_ijgk)
}
