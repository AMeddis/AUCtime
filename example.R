###example of estimating covariate-specific time dependent ROC curve and AUC

load("exdata.Rdata")
###########################
#Y:biomarker (negative bionomial)
#X:covariate (2 levels)
#status: 1 if failed, 0 if censored
#fail_time: failure time
#cluster: id of the cluster
###########################
head(exdata)
##commenges-andersen score test (to test homogeneity of failure times across clusters)
library(frailtyEM)
ca_test<-emfrail(Surv(fail_time,status)~ Y + X + cluster(cluster), data=exdata)$ca_test
ca_test
#pvalue<<0.05--> heterogeneity

#time of failure
summary(exdata$fail_time)
#time where to evaluate the AUC
timec<-quantile(exdata$fail_time,c(0.3,0.5,0.8))

##estimate the covariate-specific time dependent AUC
#without specifying teval, the AUC is estimated in all the failure times
AUCobj<-AUCtime(formula.cox="Surv(fail_time,status)~ Y + X + cluster(cluster)",formula.marker="Y~X",
                dati=exdata, status=exdata$status,marker.distr="NB",teval=timec,cov.names="X",c.values=1)
#estimated AUC
plot(AUCobj$ftimes, AUCobj$AUC,ylim=c(0.5,1),xlab="time", ylab="AUC")
#ROC curve at time=timec[2] of failure times
plot(AUCobj$FPR[2,],AUCobj$TPR[2,],type="l",lty=3,
     ylab="TPR",xlab="FPR",main="covariate-specific ROC(t=t*)")
abline(0,1)
legend(0.7,0.1,c(paste("AUC(t*)=",round(AUCobj$AUC[2],digits=3),"")),lty=3)


### conf interval for the AUC with bootstrap
#define the statistic
AUC<-function(n.data){AUCtime(formula.cox="Surv(fail_time,status) ~ Y + X + cluster(cluster)",
                              formula.marker="Y ~ X ", 
                              dati=n.data, status=n.data$status, marker.distr="NB",
                              cov.names = "X",teval=timec, c.values = 1)}

#parametric bootstrap(R=number of replications)

boot_AUC<-cp_bootAUC(data=exdata, clust.name="cluster", R=500, statistic=AUC, seed=110789)
#percentile 95% c.i.
infAUC<-sapply(1:ncol(boot_AUC$bootv),function(x){quantile(boot_AUC$bootv[,x],0.025)})
supAUC<-sapply(1:ncol(boot_AUC$bootv),function(x){quantile(boot_AUC$bootv[,x],0.975)})
tab<-cbind(boot_AUC$tv,infAUC,supAUC)
colnames(tab)<-c("AUC","AUC_l","AUC_u")
tab
