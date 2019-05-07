# Segmentation of multivariate time series with a Multiple Regression
# model with a Hidden Logistic Process (MRHLP).
#
## Please cite the following papers for this code:
#
# @article{Chamroukhi-MRHLP-2013,
# 	Author = {F. Chamroukhi and D. Trabelsi and S. Mohammed and L. Oukhellou and Y. Amirat},
# 	Journal = {Neurocomputing},
# 	Month = {November},
# 	Pages = {633--644},
# 	Publisher = {Elsevier},
# 	Title = {Joint segmentation of multivariate time series with hidden process regression for human activity recognition},
# 	Volume = {120},
# 	Year = {2013},
# 	note = {},
# 	url  = {https://chamroukhi.com/papers/chamroukhi_et_al_neucomp2013b.pdf}
# 	}
#
# @article{Chamroukhi-FDA-2018,
#  	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
#  	Author = {Faicel Chamroukhi and Hien D. Nguyen},
#  	Note = {DOI: 10.1002/widm.1298.},
#  	Volume = {},
#  	Title = {Model-Based Clustering and Classification of Functional Data},
#  	Year = {2019},
#  	Month = {to appear},
#  	url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
#     }
#
# @article{chamroukhi_et_al_NN2009,
# 	Address = {Oxford, UK, UK},
# 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
# 	Date-Added = {2014-10-22 20:08:41 +0000},
# 	Date-Modified = {2014-10-22 20:08:41 +0000},
# 	Journal = {Neural Networks},
# 	Number = {5-6},
# 	Pages = {593--602},
# 	Publisher = {Elsevier Science Ltd.},
# 	Title = {Time series modeling by a regression approach based on a latent process},
# 	Volume = {22},
# 	Year = {2009},
# 	url  = {https://chamroukhi.com/papers/Chamroukhi_Neural_Networks_2009.pdf}
# 	}
#

rm(list = ls())
source("R/FData.R")
source("R/ModelMRHLP.R")
source("R/enums.R")
source("R/ModelLearner.R")


# Building matrices for regression

# Toy multivariate time series with regime changes
# Y = cbind(c(rnorm(100, mean = 0), rnorm(120, mean = 7), rnorm(200, mean = 4), rnorm(100, mean = -1), rnorm(150, mean = 3.5)),
#            c(rnorm(100, mean = 1), rnorm(120, mean = 5), rnorm(200, mean = 6), rnorm(100, mean = -2), rnorm(150, mean = 2)),
#            c(rnorm(100, mean = -2), rnorm(120, mean = 10), rnorm(200, mean = 8), rnorm(100, mean = 0), rnorm(150, mean = 5)))
# X = matrix(seq(from = 0, to = 1, length.out = nrow(Y)), nrow = 1)

# Toy time series with regime changes
load("data/simulatedTimeSeries.RData")

# Some real time series with regime changes
# load("data/realTimeSeries.RData")

fData <- FData$new()
fData$setData(X, Y)

K <- 5 # number of regimes (mixture components)
p <- 3 # dimension of beta (order of the polynomial regressors)
q <- 1 # dimension of w (order of the logistic regression: to be set to 1 for segmentation)
variance_type <- variance_types$hetereskedastic

modelMRHLP <- ModelMRHLP(fData, K, p, q, variance_type)

# setting the model options
n_tries <- 1
max_iter = 1500
threshold <- 1e-6
verbose <- TRUE
verbose_IRLS <- FALSE


####
# EM Algorithm
####
# 1. running the em algorithm giving mixModel (data, and model itself) and the modelOptions
solution <- EM(modelMRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS)

# show the results
solution$plot()



####
# CEM Algorithm
####
#solution <- CEM(modelMRHLP, n_tries, max_iter, threshold, verbose, verbose_IRLS)
#solution$plot()
