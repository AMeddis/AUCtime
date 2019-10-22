##-----------------------------function AUCtimeloc----------------------------------------------------#
#Estimatng covariare-specific time dependent AUC for clustered survival data (clustered failure times)#
#  using shared frailty model for the hazard and a location model for the biomarker distribution      #
#-----------------------------------------------------------------------------------------------------#


#@importFrom frailtyEM emfrail, predict


#@param formula.cox:  a formula that contains on the left hand side an object of the type \code{Surv}
#and on the right hand side a \code{+cluster(id)} statement.
#@param formula.marker: regression formula for the biomarker
#@param dati:  a \code{data.frame} in which the formula argument can be evaluate
#@status: vector which defines the status of event (1 event, 0 censored)
#@param frailty_dist: specify the distribution of the frailty ("gamma","pvf","stable")
#@param teval: times where to evaluate the AUC, if NULL it is evaluated in all the failure times
#@param cov.names: the name of the covariate to adjust for 
#@param c.values: value of the covariate specific


#@return object is an object that contains this fields:
#cutoffs: vector of cut-off of the ROC curve
#ftimes: vector of time where the AUC is evaluated
#TPR: matrix of TPR(t,y) for each time(row) and cut-off(column)
#FPR: matruc of FPR (t,y) for each time(row) and cut-off(column)
#AUC: vector of the AUC at each time t
#sh.model : emfrail object with informations on the shared frailty model (see in details in the frailtyEM package)
#mark_coef: coefficients for the biomarker with coeff of the regression 
#survP: marginal survival function S(t)
#survM: marginal survival function respect to the frailty term S(t |Y,X)

#@estimation
#shared frailty model: the partial likelihood and the EM algorithm is used.
#The algorithm is detailed in the package frailtyEM.
#location model as detailed by Song ans Zhou


AUCtimeloc<-function( formula.cox, formula.marker , dati, status, frailty_dist="gamma", teval=NULL,
                      cov.names=NULL, c.values=NULL )
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
  
  
  #-------------for the biomarker
  #linear regression
  markmod<-lm(formula.marker, data=dati)
  beta<-coef(markmod)
  #calculate the residuals centered in the specific cov value
  cval<-1
  res<-dati$marker - beta[1]
  
  if(!is.null(cov.names))
  { cval<- c(1,c.values)
  res<- res - beta[2]*(dati[ ,cov.names] - c.values)}
  
  dimv<-length(res)
  #dimv<-nrow(res)
  
  #parameters 
  beta<-c(beta)
  
  
  ###fit frailty model
  model<-frailtyEM::emfrail(data=dati, formula=formula.cox, 
                            emfrail_dist(dist=frailty_dist), emfrail_control(ca_test=FALSE,lik_ci = FALSE))
  ##time of event
  times<-predict(model,quantity = "survival", type = "marginal",lp=c(0),conf_int=NULL)$time
  
  
  ##compute marginal (respect to frailty) Survival probabilities for each value of marker
  if(!is.null(cov.names))
  { 
    #create data frame where estimate the surv function
    #marker.matrix <- matrix(res[,1],nrow=dimv,ncol=1)
    marker.matrix <- matrix(res,nrow=dimv,ncol=1)
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
      marker.matrix <- matrix(res[j,1],nrow=1,ncol=1)
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
  
  
  #Compute the TPR and FPR with 
  #TPR(c,t)=P(CTC >= c|X,T<=t) 
  #FPR(c,t)=P(CTC >= c|X,T>t) 
  
  #numerator of TPR and FPR
  #we have to sum over the observations 
  #where res>obsvalue (osbvalue is c in the formula above)
  num<-do.call(cbind,lapply(obsvalue,FUN=function(i){
    k<-which(res >= i)
    numT<-apply(cbind((1-surv[,k]),rep(0,dimt)),1,sum)
    numF<-apply(cbind(surv[,k],rep(0,dimt)),1,sum)
    c(numT,numF)}))
  
  #den we sum over all the values
  denT<-apply((1-surv),1,sum)
  denF<-apply(surv,1,sum)
 
  
  
  #TPR and FPR (we add 1 and 0 to the TPR and FPR)
  TPR<-cbind(rep(1,dimt),num[1:dimt,]/rep(denT,dimobsv),rep(0, dimt))
  FPR<-cbind(rep(1,dimt),num[(dimt+1):(2*dimt),]/rep(denF,dimobsv),rep(0, dimt))
  
  
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

