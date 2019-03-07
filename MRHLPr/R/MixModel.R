source("R/dataset.R")
source("R/utils.R")

MixModel <- setRefClass(
  "MixModel",
  contains = "MyData",
  # Define the fields
  fields = list(K="numeric", # number of regimes (mixture components)
                p="numeric", # dimension of beta' (order of the polynomial regressors)
                q="numeric" # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
                ),
  methods = list(

  )
)

MixModel<-function(mixData, K,p,q){
  new("MixModel",x=mixData$x, y=mixData$y, m=mixData$m, d=mixData$d, K=K, p=p, q=q)
}


