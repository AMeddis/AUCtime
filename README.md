# AUCtime
Function useful to estimate covariate-specific time dependent ROC curve and AUC for a candidate biomarker in case of correlated censored failure times.

We consider a shared frailty model for the hazard and a parametric model for the biomarker given the covariate. In particular, we assume a negative binomial distribution
for the biomarker and an umbalanced covariate with 2 levels, as in our motivating example (CTCs count in non metastatic breast cancer).

We provide:

- exdata: data set of 5 variables. Y is the biomarker, X the covariate, fail_time the times of failure, cluster the cluster id for each subject
 and status is 1 if the subject occurs the event, 0 otherwise (censored data).
 
- AUCtime: function to estimate the covariate-specific time dependent ROC curve and AUC specifying the formula of the shared frailty model 
 and the biomarker model. It returns the covariate-specific time dependent TPR and FPR for any treshold, the AUC(t*) if t* is specified, otherwise 
 the AUC at all the observed failure times is provided.
 
- parbootNICS: functions used for the confidence interval obtained by parametric bootstrap

  - parboot: function to generate the bootstrap data set from the estimated parameters of the original data set
  
  - cp_bootAUC: function to obtain the estimated AUC for each boostrap dataset
 
 - example: an example of how to estimate the covariate-specific time dependent AUC with its confidence interval. The exdata database is 
 used for this example.


