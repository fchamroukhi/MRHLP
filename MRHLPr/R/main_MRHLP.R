"
% Segmentation of multivariate time series with a Multiple Regression
% model with a Hidden Logistic Process (MRHLP).
%
%% Please cite the following papers for this code:
%
% @article{Chamroukhi-MRHLP-2013,
% 	Author = {F. Chamroukhi and D. Trabelsi and S. Mohammed and L. Oukhellou and Y. Amirat},
% 	Journal = {Neurocomputing},
% 	Month = {November},
% 	Pages = {633--644},
% 	Publisher = {Elsevier},
% 	Title = {Joint segmentation of multivariate time series with hidden process regression for human activity recognition},
% 	Volume = {120},
% 	Year = {2013},
% 	note = {},
% 	url  = {https://chamroukhi.com/papers/chamroukhi_et_al_neucomp2013b.pdf}
% 	}
%
% @article{Chamroukhi-FDA-2018,
%  	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
%  	Author = {Faicel Chamroukhi and Hien D. Nguyen},
%  	Note = {DOI: 10.1002/widm.1298.},
%  	Volume = {},
%  	Title = {Model-Based Clustering and Classification of Functional Data},
%  	Year = {2019},
%  	Month = {to appear},
%  	url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
%     }
%
% @article{chamroukhi_et_al_NN2009,
% 	Address = {Oxford, UK, UK},
% 	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
% 	Date-Added = {2014-10-22 20:08:41 +0000},
% 	Date-Modified = {2014-10-22 20:08:41 +0000},
% 	Journal = {Neural Networks},
% 	Number = {5-6},
% 	Pages = {593--602},
% 	Publisher = {Elsevier Science Ltd.},
% 	Title = {Time series modeling by a regression approach based on a latent process},
% 	Volume = {22},
% 	Year = {2009},
% 	url  = {https://chamroukhi.com/papers/Chamroukhi_Neural_Networks_2009.pdf}
% 	}
%
"
rm(list = ls())
source("R/dataset.R")
source("R/enums.R")
source("R/MixModel.R")
source("R/ModelOptions.R")
source("R/ModelLearner.R")


# loading and setting the data
fileName = "R/datasets/simulated_time_series.mat"
mixData <- MyData$new()
mixData$setDataFromMat(fileName)

# setting the model
K <- 5; # number of regimes (mixture components)
p <- 3; # dimension of beta' (order of the polynomial regressors)
q <- 1; # dimension of w (ordre of the logistic regression: to be set to 1 for segmentation)
mixModel <- MixModel(mixData,K,p,q)

# setting the model options
n_tries <- 1 # number of tries EM/CEM to run
max_iter <- 1500 # maximum number of iteration in the EM/CEM algorithm
threshold <- 1e-6 # threshold to check the concergence
verbose <- TRUE # verbose the EM/CEM algorithm
verbose_IRLS <- FALSE # verbose the IRLS algorithm
modelOptions <- ModelOptions(n_tries, max_iter, threshold, verbose, verbose_IRLS, variance_types$hetereskedastic)

####
# EM Algorithm
####
# 1. running the em algorithm giving mixModel (data, and model itself) and the modelOptions
solution <- EM(mixModel, modelOptions)
# 2. getting the mixture parameters solution
mixParamSolution <- solution[[1]]
# 3. getting the mixture stats solution
mixStatsSolution <- solution[[2]]

# show the results
mixStatsSolution$showDataClusterSegmentation(mixModel, mixParamSolution)
####
# CEM Algorithm
####
#solution <- CEM(mixModel, modelOptions)
#mixParamSolution <- solution[[1]]
#mixStatsSolution <- solution[[2]]
#mixStatsSolution$showDataClusterSegmentation(mixModel, mixParamSolution)

