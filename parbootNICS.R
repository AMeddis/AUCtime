#function to provide the parametric bootstrap for clustered data with NICS
##parboot function to generate da bootstrap data set 
#cp_bootAUC function to estimate the AUC in each bootstrap data set

###1.parboot function
#@param data: dataset for the resampling
#@param clust.name:  name of the variable that identifies the cluster
#@param est.mark_coef: estimate coeff of the biomarker (regression coeff and parameters)
#@param mod: estimated shared frailty model
#@param Cens: object obtained from cp_bootAUC that is the Kaplan-Meier estimator of the censoring distribution

#@importFrom dplyr group_by, count
#@importFrom frailtyEM predict.emfrail

#@return object is a dataframe with:
#"cov.names":covariate value
#status: 1 if failed, 0 if censored
#"mar.names": marker value
#"time.names": vector of failure times
#"clust.name": vector of the cluster id

parboot<-function(data,clust.name, est.mark_coef, mod, Cens)
{
  library(dplyr)
  library(frailtyEM)
  #save the name of the variables from the formula
  cov.names<-as.character(attr(terms.formula(mod$formula), "variables")[4])
  mar.names<-as.character(attr(terms.formula(mod$formula), "variables")[3])
  time.names<-as.character(attr(terms.formula(mod$formula), "variables")[[2]][2])
  stat.names<-as.character(attr(terms.formula(mod$formula), "variables")[[2]][3])
  data$family<-eval(parse(text=paste("data$",clust.name,"", sep="")))
  
  ##1. sample the cluster sizes with replecement
  # number of cluster
  n_cluster<-length(table(data$family))
  #save the clusters sizes
  size_cluster<-(data %>% group_by(family) %>% count(n()))$n
  #sample with replecement
  clusterss<-sample(size_cluster, replace=TRUE)
  #new cluster sample sizes
  n_new<-sum(clusterss)
  
  
  #generate U from the gamma distribution with the estmated teta
  est.teta<-exp(mod$logtheta)
  new_Uc<-rgamma(n_cluster,est.teta, est.teta)
  
  #resample n_new covariate Z from original data
  indZ<-sample(c(1:sum(size_cluster)),n_new, replace=TRUE)
  Z<-data[indZ,cov.names]
  
  #we assume a negative binomial biomarker
  #parameters for neg binomial distr
  media<-function(z){exp(t(est.mark_coef[1:2]) %*%c(1,z))}
  d<-est.mark_coef[4]
  
  #generate Y given Z
  cov_val<-as.numeric(names(table(Z)))
  Y<-rep(0,n_new)
  for(j in 1:length(cov_val))
  {Y[which(Z==cov_val[j])]<-rnbinom(length(which(Z==cov_val[j])) ,mu=(media(cov_val[j])),size=d)}
  
  ##generate the time of events 
  #estimate coef of the shared frailty model
  est.Cox_coef<-mod$coefficients 
  #newdata for the prediction of Lambda0
  dframe<-function(y,k,z){
    datamatrix <- matrix(c(y,z),nrow=clusterss[k],ncol=2)
    colnames(datamatrix) <- c(mar.names,cov.names)
    data.frame(datamatrix )}
  
  #function to generate time of event in each cluster i
  a<-function(i){
    #z and y for the i-th cluster  
    ind_cl<-ifelse(i==1,1,sum(clusterss[1:(i-1)])+1)
    z<-Z[ind_cl:(ind_cl + (clusterss[i]-1))]
    y<-Y[ind_cl:(ind_cl + (clusterss[i]-1))]
    
    #the conditional cumhazard function
    Lambda<-predict(mod,quantity = "cumhaz", type = "conditional", newdata = dframe(y,i,z),conf_int=NULL)
    
    #generate times-of-event inverting the the estimated cumulative distribution function of T|Y,X,U
    V<-runif(clusterss[i])
    tempi<-as.numeric(sapply(1:clusterss[i],function(x){ 
      tmpt<-c(Lambda[[x]]$time,max(Lambda[[x]]$time))
      cdf=1-exp(-Lambda[[x]]$cumhaz*new_Uc[i])  
      cdf<-c(cdf,1)
      tmpt[min(which(cdf >=V[x]))]}))
    #create the data.frame for cluster i
    data.frame("Z"=z, "Y"=y, "time"=tempi, "family"=rep(i, clusterss[i]))
    }
  
  ##boostrap data
  new.data<-do.call(rbind,lapply(seq(1:n_cluster),FUN=a))
 
  ##Ad censoring
  #surv funct for censoring time
  censtime<- eval(parse(text=paste("unique(sort(data$",time.names,"[data$",stat.names,"==0]))")))
  survCens <- predict(Cens,times=censtime)
  #it is piecewise function
  Jumps <- -diff(c(1,survCens,0))
  #censoring time obtained sampling from the censoring time of the original data set 
  #with probability obtained by the estimated distribution of the censoring
  ctime <- sample(x=c(censtime,max(censtime)+1),size=n_new,replace=TRUE,prob=Jumps)
  #definition of the status
  new.data$status<-as.numeric(new.data$time <= ctime)
  new.data$time <- pmin(new.data$time, ctime)
  
  #bootstrap data set where to evaluate the AUC
  nd<-eval(parse(text=paste("data.frame(",cov.names,"=new.data$Z,", stat.names,"=new.data$status,",
                            mar.names,"= new.data$Y,", time.names,"=new.data$time,", 
                            clust.name,"=new.data$family)")))
  
  nd
  
}  



###2.cp_bootAUC function 
#@param data: original data set where to apply the bootstrap
#@param clust.name: name of the variable that identifies the cluster
#@param R: number of replication of bootstrap
#@param statistic: function that provide the statistic we want to calculate the bootstrap for
#@param seed: seed to have reproducible results

#@importFrom prodlim prodlim

#@return object is a list with
#tv: is the estimated AUC for the original data set
#bootv: a matrix with the estimated AUC for each bootstrap data sets (R rows * dimt columns)
#parboot: data frame of the estimated parameters from the boostrap data sets
#partv: data frame of the estimated parameters from the original data set

cp_bootAUC<-function(data, clust.name, R, statistic, seed=1805)
{
  library(prodlim)
  #true value
  tmp<-statistic(data)
  t0<-tmp$AUC
  #we need the estimated shared frailty model and the parameters for the biomarker
  est.mark_coef<-tmp$mark_coef
  mod<-tmp$sh.mod
  
  #estimated censoring distribution
  f<-as.character(attr(terms.formula(mod$formula),"variables")[2])
  formula<-eval(parse(text=paste("",f,"~1", sep="")))
  Cens<- prodlim(formula , data, reverse=TRUE)
  
  set.seed(seed)
  
  #create R bootstrap data set where to calculate the statistic
  #we obtain AUC for each bootstrap data set and the respective estimated parameters
  bt<- do.call(cbind,lapply(seq(1:R), FUN=function(i){
    newdata<-parboot(data, clust.name, est.mark_coef, mod, Cens)
    ris<-statistic(newdata)
    list("AUC"=ris$AUC,"theta"=exp(ris$sh.mod$logtheta), "beta0"=ris$sh.mod$coefficients[1],
         "beta1"=ris$sh.mod$coefficients[2], "mu"=ris$mark_coef[3], "d"=ris$mark_coef[4])
    
  }))

  #create the data frame with the estimated parameters of the bootstrap data sets
  paramboot<-do.call(rbind, lapply(1:R,FUN=function(k){data.frame("beta0"=bt[,k]$beta0,"beta1"=bt[,k]$beta1,
                                                                "theta"=bt[,k]$theta,"mu"=bt[,k]$mu, "d"=bt[,k]$d)}))
  #create the data frame with the estimated parameter of the true value
  partv<-data.frame("beta0"=mod$coefficients[1],"beta1"=mod$coefficients[2],
                    "theta"=exp(mod$logtheta),"mu"=est.mark_coef[3], "d"=est.mark_coef[4])
  #create the data frame with the bootstrap AUC values
  AUC<-do.call(rbind, lapply(1:R,FUN=function(k){bt[,k]$AUC}))
  
  #output
  list("tv"=t0, "bootv"=AUC, "parboot"=paramboot, "partv"=partv)
  
}
