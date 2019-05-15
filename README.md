# R codes for Segmentation of multivariate time series with a Multiple Regression model with a Hidden Logistic Process (MRHLP).: 

<< Note that the codes are also provided in Matlab and Python >>

R codes written by

**Faicel Chamrouckhi**
&
**Marius Bartcus**
&
**Florian Lecocq**

firstname.lastname@unicaen.fr


When using this code please cite the following papers : The two first ones concern the model and its use in clusterng and the last ones concern the model and its use in discrimination.

The script to run is: main_MixFRHLP_EM.R. The main function is in ModelLearner.R:
EM(...) for the EM algorithm

```
 @article{Chamroukhi-MRHLP-2013,
	Author = {F. Chamroukhi and D. Trabelsi and S. Mohammed and L. Oukhellou and Y. Amirat},
	Journal = {Neurocomputing},
	Month = {November},
	Pages = {633--644},
	Publisher = {Elsevier},
	Title = {Joint segmentation of multivariate time series with hidden process regression for human activity recognition},
	Volume = {120},
	Year = {2013},
	note = {},
	url  = {https://chamroukhi.com/papers/chamroukhi_et_al_neucomp2013b.pdf}
	}

@article{Chamroukhi-FDA-2018,
 	Journal = {Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery},
 	Author = {Faicel Chamroukhi and Hien D. Nguyen},
 	Note = {DOI: 10.1002/widm.1298.},
 	Volume = {},
 	Title = {Model-Based Clustering and Classification of Functional Data},
 	Year = {2019},
 	Month = {to appear},
 	url =  {https://chamroukhi.com/papers/MBCC-FDA.pdf}
    }

@article{chamroukhi_et_al_NN2009,
	Address = {Oxford, UK, UK},
	Author = {Chamroukhi, F. and Sam\'{e}, A. and Govaert, G. and Aknin, P.},
	Date-Added = {2014-10-22 20:08:41 +0000},
	Date-Modified = {2014-10-22 20:08:41 +0000},
	Journal = {Neural Networks},
	Number = {5-6},
	Pages = {593--602},
	Publisher = {Elsevier Science Ltd.},
	Title = {Time series modeling by a regression approach based on a latent process},
	Volume = {22},
	Year = {2009},
	url  = {https://chamroukhi.com/papers/Chamroukhi_Neural_Networks_2009.pdf}
	}
```


### THE SHORT DESCRIPTION OF EACH R FILE. For more detailed description, please see the individual files

1) ModelLearner _Contains the two functions of the EM and the CEM algorithm._
2) FData _Contains the dataset object._                        
3) ModelMRHLP _The ModelMixRHLP class containts the data object and the model settings (number of clusters, the number of regimes, the degree of polynomials, the order of the logistic regression)_
4) ParamMRHLP _Initializes and updates (the M-step) the model parameters (parameters of the logistic process for each of the clusters, polynomial regression coefficients, and the variances for the regmies for each cluster)._
5) StatMRHLP _Calculates mainly the posterior memberships (E-step), the loglikelihood, the parition, different information criterias BIC, ICL, etc_
6) enums _Used to enumerate the variance type (heteroskedastic or homoscedastic)_
7) utils _Contains some helping functions_
8) myKmeans _Kmeans algorithm_
9) FittedMRHLP _The learned model_


