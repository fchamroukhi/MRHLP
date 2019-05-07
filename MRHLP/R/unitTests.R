rm(list = ls())
source("R/dataset.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/enums.R")
source("R/RegressionDesigner.R")
source("R/MixParam.R")

testDataSets <- function(){
  fileName = "R/datasets/simulated_time_series.mat"
  mixData <- MyData$new()
  mixData$setDataFromMat(fileName)
  return(mixData)
}

mixData <- testDataSets()


testMixModel <- function(){
  fileName = "R/datasets/simulated_time_series.mat"
  mixData <- MyData$new()
  mixData$setDataFromMat(fileName)

  # setting the model
  K <- 5; # number of regimes (mixture components)
  p <- 3; # dimension of beta' (order of the polynomial regressors)
  q <- 1; # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
  mixModel <- MixModel(mixData,K,p,q)

  return(mixModel)
}

mixModel <- testMixModel()

testMixOptions <- function(){
  # setting the model options
  n_tries <- 1 # number of tries EM/CEM to run
  max_iter <- 1500 # maximum number of iteration in the EM/CEM algorithm
  threshold <- 1e-6 # threshold to check the concergence
  verbose <- TRUE # verbose the EM/CEM algorithm
  verbose_IRLS <- FALSE # verbose the IRLS algorithm
  modelOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$hetereskedastic)
  return(modelOptions)
}

modelOptions <- testMixOptions()

testRegressionDesigner <- function(){
  mixModel <- testMixModel()
  phi <- RegressionDesigner$new()
  phi$setPhi1(mixModel$x,mixModel$p,mixModel$q)
  return(phi)
}

regDesigner <- testRegressionDesigner()





testInitializeParameters <- function(){
  fileName = "R/datasets/simulated_time_series.mat"
  mixData <- MyData$new()
  mixData$setDataFromMat(fileName)

  # setting the model
  K <- 5; # number of regimes (mixture components)
  p <- 3; # dimension of beta' (order of the polynomial regressors)
  q <- 1; # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
  mixModel <- MixModel(mixData,K,p,q)

  try_EM <- 1
  phi <- RegressionDesigner$new()
  phi$setPhi1(mixModel$x,mixModel$p,mixModel$q)

  # Initialization
  mixParam <- MixParam(mixModel, modelOptions)
  mixParam$initRegressionParam(phi$XBeta, mixModel$y, mixModel$K, variance_types$hetereskedastic, try_EM)
  return(mixParam)
}

mixParam <- testInitializeParameters()
