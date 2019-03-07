rm(list = ls())
source("R/dataset.R")

testDataSets <- function(){
  fileName = "R/datasets/simulated_time_series.mat"
  mixData <- MyData$new()
  mixData$setDataFromMat(fileName)
  return(mixData)
}

mixData = testDataSets()
