source("R/enums.R")
source("R/utils.R")
MixParam <- setRefClass(
  "MixParam",
  fields = list(
    Wg = "matrix", # Wg = (Wg1,...,w_gK-1) parameters of the logistic process:
                  # matrix of dimension [(q+1)x(K-1)] with q the order of logistic regression.
    betag = "array", # betag = (beta_g1,...,beta_gK) polynomial regression coefficient vectors: matrix of
                     # dimension [(p+1)xK] p being the polynomial  degree.
    sigmag = "array" # sigma_g = (sigma_g1,...,sigma_gK) : the variances for the K regmies. vector of dimension [Kx1]
  ),
  methods = list(
    initRegressionParam = function(XBeta, y, K, variance_type, try_algo){
      "
        Initializes the regressions parameters for the MRHLP model: the
        regressions coefficients vector and the variance, for each component.
      "
       n <- nrow(XBeta) # m
       m <- ncol(XBeta) # P

       if (try_algo==1){ #  uniform segmentation into K contiguous segments, and then a regression
          zi <- round(n/K)-1

          s <- 0

          for (k in 1:K){
            i <- (k-1)*zi+1
            j <- k*zi

            yk <- y[i:j,]
            Xk <- XBeta[i:j,]

            muk <- Xk %*% betag[,,k]
            sk <- t(yk - muk) %*% (yk - muk)

            if (variance_type == variance_types$homoskedastic){
              s <- s + sk
              sigmag <<- s/n
            }
            else{
              sigmag[, ,k] <<- sk / length(yk)
            }
          }
       }
       else{ # random segmentation into contiguous segments, and then a regression
         Lmin <- m+1 #round(m/K) #nbr pts min into one segment
         tk_init <- zeros(K,1)
         K_1 <- K
         for (k in 2:K) {
           K_1 <- K_1-1;
           #temp <- (tk_init[k-1] + Lmin) : (m - (K_1*Lmin))

           temp <- tk_init[k-1] + Lmin : (n - (K_1*Lmin) - tk_init[k-1])

           ind <- sample(length(temp));
           tk_init[k] <- temp[ind[1]]
         }
         tk_init[K+1] <- n

         s <- 0
         for (k in 1:K){
           i <- tk_init[k] + 1
           j <- tk_init[k+1]

           yk <- y[i:j,]
           Xk <- XBeta[i:j,]


           betag[,,k] <<- solve(t(Xk)%*%Xk)%*%t(Xk)%*%yk

           muk <- Xk %*% betag[,,k]
           sk <- t(yk - muk) %*% (yk - muk)

           if (variance_type == variance_types$homoskedastic){
             s <- s + sk
             sigmag <<- s/n
           }
           else{
             sigmag[, ,k] <<- sk / length(yk)
           }
         }
       }
    }

  )
)

MixParam<-function(mixModel, options){
  #mixModel <- mixModel
  Wg <- rand(mixModel$q+1, mixModel$K-1)
  betag <- array(NA, dim=c(mixModel$p+1, mixModel$d, mixModel$K))
  if (options$variance_type == variance_types$homoskedastic){
    sigmag <- matrix(NA, mixModel$d, mixModel$d)
  }
  else{
    sigmag <- array(NA, dim = c(mixModel$d, mixModel$d, mixModel$K))
  }
  new("MixParam", Wg=Wg, betag=betag, sigmag=sigmag)
}
