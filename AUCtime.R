#######function AUCtime######
#estimatng covariare-specific time dependent AUC for clustered survival data (clustered failure times)
#using frailty model for the hazard and a parametric model for the biomarker distribution for both
#discrete biomarker with negative bionomial distribution and continuous biomarker (normal distribution)


#@importFrom frailtyEM emfrail, predict
#@importFrom pscl glm.nb

#@param formula.cox:  a formula that contains on the left hand side an object of the type \code{Surv}
#and on the right hand side a \code{+cluster(id)} statement.
#@param formula.marker: regression formula for the biomarker
#@param dati:  a \code{data.frame} in which the formula argument can be evaluate
#@status: vector which defines the status of event (1 event, 0 censored)
#@param frailty_dist: specify the distribution of the frailty ("gamma","pvf","stable")
#@param teval: times where to evaluate the AUC, if NULL it is evaluated in all the failure times
#@param marker.distr: it defines the model for the biomarker, 
#if "NB" discrete biomarker with negative binomial distribution, otherwise continuous biomarker (normal distribution)
#@param lowerP.bmk: defines the stop rules for the sum over the biomarker value in the formula of TPR and FPR.
#if the P(Y=y)<lowerP.bmk than STOP: the contribute to the sum is negligible
#@param cov.names: the name of the covariate to adjust for 
#@param c.values: value of the covariate specific

#@return object is an object that contains this fields:
#cutoffs: vector of cut-off of the ROC curve
#ftimes: vector of time where the AUC is evaluated
#TPR: matrix of TPR(t,y) for each time(row) and cut-off(column)
#FPR: matruc of FPR (t,y) for each time(row) and cut-off(column)
#AUC: vector of the AUC at each time t
#sh.model : emfrail object with informations on the shared frailty model (see in details in the frailtyEM package)
#mark_coef: coefficients for the biomarker with coeff of the regression and the distribution parameters:
#discrete biomarker: negativ binomial with mean and dispersion parameter
#continuous biomarker: normal distribution with mean and variance
#survP: marginal survival function S(t)
#survM: marginal survival function respect to the frailty term S(t |Y,X)

#@estimation
#shared frailty model: the partial likelihood and the EM algorithm is used.
#The algorithm is detailed in the package frailtyEM.
#parametric model: classical regression model are used

#@confidence interval
#a parametric bootstrap is used for 95% percentile confidence interval


AUCtime<-function( formula.cox, formula.marker , dati, status, frailty_dist="gamma", teval=NULL,marker.distr="NULL", 
                   lowerP.bmk= 10e-6 ,cov.names=NULL, c.values=NULL )
{
  library(frailtyEM)
  library(dplyr)
  
  
  if(is.character(formula.cox)) {
    formula.cox <- as.formula(formula.cox)
  }
  if(is.character(formula.marker)) {
    formula.marker <- as.formula(formula.marker)
  }
  
  
  marker <- as.character(attr(terms.formula(formula.marker), "variables")[2])
  dati$marker<-dati[ , c(marker)]
  dati$status<-status
  
  
  if(!is.null(cov.names))
  {if(is.null(c.values) || length(c.values) != length(cov.names))
    stop('no covariate value specified')
  }
  
  
  #save all the observed values of the biomarker
  #that define the cut-off for TPR and FPR
  obsvalue <- as.numeric(names(table(dati$marker)))
  dimobsv<-length(obsvalue)
  max_value<-max(obsvalue)
  
  
  
  if(marker.distr=="NB")
  {
    library(pscl)
    #discrete biomarker Y (negative bionomial distribution)
    
    #for the TPR and FPR we want to approximate \int_{y}^{\infty} --> sum_{0}^{M}
    #M=maxvalue+100, (it is negative binomial, namely zero inflated)
    value<-seq(0,max_value+100,by=1) 
    
    #negative binomial regression
    NB<-glm.nb(formula.marker,data=dati)
    beta<-coef(NB)
     
    cval<-1
    if(!is.null(cov.names))
    { cval<- c(1,c.values)}
    #mean 
    mu<-exp(t(beta) %*% cval)
    
    #P(Y=y|X) (dispersion parameter= NB$theta)
    P<-dnbinom(value,mu=mu,size=NB$theta)
    #STOP criteria for the sum for TPR and FPR
    excl<-which(P<lowerP.bmk & value>max_value)
    #save in value the biomarker values for the sum for TPR and FPR
    if(length(excl)!=0)
    {value<-value[-excl]
    P<-P[-excl]
    }
   
    
    #in beta we save the coeff of negative binomial
    #beta: coeff of the neg binomial regression
    #mu, NB$theta: parameters for the negative binomial distribution
    #(these parameters are needed for the parametric bootstrap)
    beta<-c(beta,mu,NB$theta)
  }
  else{
    ##continuous biomarker Y (normal distribution)
  
    #linear regression
    markmod<-lm(formula.marker, data=dati)
    beta<-coef(markmod)
    
    cval<-1
    if(!is.null(cov.names))
    { cval<- c(1,c.values)}
    
    #distribution parameters  
    mu<-(t(beta) %*% cval)
    sigma<-summary(markmod)$sigma
    
    #for the TPR and FPR we want to approximate \int_{y}^{\infty} --> sum_{m}^{M}
    #we generate a normal distribution with the estimated parameters 
    #with sample sizes=2000 (if the sample size is smaller)
    ss<-ifelse(nrow(dati)<2000,2000,nrow(dati))
    value<-seq(min(obsvalue)-2*sigma,max_value+2*sigma, length=ss)
    
    #P(Y=y | X)
    P<-dnorm(value,mu,sigma)
    #STOP criteria for the sum for TPR and FPR
    excl<-which(P<lowerP.bmk)
    #save in value the biomarker values for the sum for TPR and FPR
    if(length(excl)!=0)
    {value<-value[-excl]
    P<-P[-excl]
    }
    #in beta we save the parameters
    #beta: regression param
    #mu sigma: distribution parameters
    beta<-c(beta,mu,sigma)
    }
  
  dimv<-length(value)
  
  ###fit frailty model
  model<-frailtyEM::emfrail(data=dati, formula=formula.cox, 
                 emfrail_dist(dist=frailty_dist), emfrail_control(ca_test=FALSE,lik_ci = FALSE))
  ##time of event
  times<-predict(model,quantity = "survival", type = "marginal",lp=c(0),conf_int=NULL)$time
  
  
  ##compute marginal (respect to frailty) Survival probabilities for each value of marker

  if(!is.null(cov.names))
  { 
    #create data frame where estimate the surv function
    marker.matrix <- matrix(value,nrow=dimv,ncol=1)
    colnames(marker.matrix) <- marker
    cov.matrix<-matrix(c.values,nrow=1,ncol=length(c.values))
    colnames(cov.matrix)<-cov.names
    df<-data.frame(marker.matrix,cov.matrix)
    #survival function S(t|Y,X)
    surv<-do.call(cbind,lapply(seq(1:dimv), FUN=function(j){predict(model, quantity = "survival", type = "marginal", 
                                                                    newdata = df[j,],conf_int=NULL )$survival_m}))
    }
  else{
    #survival function S(t|Y)
    surv<-do.call(cbind,lapply(seq(1:dimv), FUN=function(j){
                                                     marker.matrix <- matrix(value[j],nrow=1,ncol=1)
                                                     colnames(marker.matrix) <- marker
                                                predict(model, quantity = "survival", type = "marginal", 
                                                     newdata = data.frame(marker.matrix),conf_int=NULL )$survival_m}))}
  
  ##time where to evaluate AUC(t)
  if(!is.null(teval))
  {
    ind<-sapply(teval,FUN=function(i){max(which(times<=i))})
    times<-teval
    surv<-surv[ind,]
  }
  dimt<-length(times)
  
  
  #####COmpute the TPR and FPR with  TPR(c,t)=P(CTC >= c | T<=t)
  
  #numerator of TPR and FPR
  num<-do.call(cbind,lapply(obsvalue,FUN=function(i){
                                         k<-min(which(value>=i))
                                         numT<-(1-surv[,k:dimv])%*%as.matrix(P[k:dimv])
                                         numF<-surv[,k:dimv]%*%as.matrix(P[k:dimv])
                                         rbind(numT,numF)}))
  
  #denominator of TPR and FPR
    denT<-(1-surv)%*%P
    denF<-surv%*%P 
  
  #TPR and FPR
   TPR<-cbind(num[1:dimt,]/rep(denT,dimobsv),rep(0, dimt))
   FPR<-cbind(num[(dimt+1):(2*dimt),]/rep(denF,dimobsv),rep(0, dimt))
  

  #function to compute AUC by trapezoidal rule
  AireTrap<-function(Abs,   # values of x-axis
                     Ord){  # values of y-axis
    nobs <- length(Abs)
    dAbs <- -(Abs[-1]-Abs[-nobs])
    mil <- (Ord[-nobs]+Ord[-1])/2
    area <- sum(dAbs*mil)
    return(area)
  }
  
  
  AUC<-sapply(1:dimt, FUN = function(i){ AireTrap(FPR[i,],TPR[i,])})
  
  #marginal survival function S(t)
  survP<-predict(model, quantity = "survival", type = "marginal", lp=c(0), conf_int = NULL)$survival_m
  if(!is.null(teval))
  {survP<-survP[ind]}
  
  ris <- list(cutoffs=c(-Inf,obsvalue,Inf), ftimes=times, TPR = TPR, FPR = FPR,
              AUC = AUC, sh.mod= model , mark_coef= beta ,survP=survP,survM=surv)
  ris
}
