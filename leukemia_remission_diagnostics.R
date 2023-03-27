## Leukemia remission data
## Diagnostic graphs and tests for PH assumption
library(survival)
library(rms)
#library(survminer) # required for ggcoxzph()
## Read in leukemia remission data
leuk=read.csv("C:/Users/andre/Desktop/Cornell/Junior Spring/STSCI4270_Survival_Analysis/Code Data/leukemia_remission_with_sex.csv")
leuk$group=as.factor(leuk$group)
leuk$sex=as.factor(leuk$sex)
## KM analysis for Rx groups
kmfit=npsurv(Surv(remtime,status)~group,data=leuk) # survfit is depreciated
## ln[-lnS] plot: treatment versus placebo
plot(kmfit,lty=c("solid","dashed"),col=c("black","blue"),
     fun="cloglog",main="ln[-ln S(t)] plot")
legend("topleft",c("treatment","placebo"),
       lty=c("solid","dashed"),col=c("black","blue"))
survplot(kmfit,loglog=TRUE)
## ln[-lnS] plot: by logWBC category
leuk$WBCgroup=cut(leuk$logwbc,breaks=c(1.4,2.4,3.1,5.1),
                  labels=c("low","medium","high"))
kmfit=npsurv(Surv(remtime,status)~WBCgroup,data=leuk)
plot(kmfit,col=c("green","black","blue"),
     fun="cloglog",main="ln[-ln S(t)] plot")
legend("topleft",c("low","medium","high"),lty=rep("solid",3),
       col=c("green","black","blue"))
survplot(kmfit,loglog=TRUE)
## ln[-lnS] plot: men versus women
kmfit=npsurv(Surv(remtime,status)~sex,data=leuk)
plot(kmfit,lty=c("solid","dashed"),col=c("black","blue"),
     fun="cloglog",main="ln[-ln S(t)] plot")
legend("topleft",c("women","men"),
       lty=c("solid","dashed"),col=c("black","blue"))
survplot(kmfit,loglog=TRUE,logt=TRUE)
##
## treatment versus placebo adjusted for logWBC
coxfit=npsurv(coxph(Surv(remtime,status)~logwbc+strata(group),data=leuk))
survplot(coxfit,loglog=TRUE,conf="none")
##
## WBC group adjusted for treatment
coxfit=npsurv(coxph(Surv(remtime,status)~group+strata(WBCgroup),data=leuk))
survplot(coxfit,loglog=TRUE,conf="none")
##
## sex adjusted for treatment and logWBC
coxfit=npsurv(coxph(Surv(remtime,status)~group+logwbc+strata(sex),data=leuk))
survplot(coxfit,loglog=TRUE,conf="none",logt=TRUE)
##
## Observed versus expected plot for group
coxfit=coxph(Surv(remtime,status)~group,data=leuk)
kmfit=npsurv(Surv(remtime,status)~group,data=leuk)
survplot(kmfit,conf.int=FALSE,lty=1,lwd=2,col="black")
newdat=data.frame(group=c("0","1"))
lines(npsurv(coxfit,newdata=newdat),col="red",lty=1,lwd=2)
##
## Observed versus expected plot for logWBC
coxfit=coxph(Surv(remtime,status)~WBCgroup,data=leuk)
kmfit=npsurv(Surv(remtime,status)~WBCgroup,data=leuk)
survplot(kmfit,conf.int=FALSE,lty=1,lwd=2,col="black")
newdat=data.frame(WBCgroup=c("low","medium","high"))
lines(npsurv(coxfit,newdata=newdat),col="red",lty=1,lwd=2)
##
## Observed versus expected plot for sex
coxfit=coxph(Surv(remtime,status)~sex,data=leuk)
kmfit=npsurv(Surv(remtime,status)~sex,data=leuk)
survplot(kmfit,conf.int=FALSE,lty=1,lwd=2,col="black")
newdat=data.frame(sex=c("0","1"))
lines(npsurv(coxfit,newdata=newdat),col="red",lty=1,lwd=2)
##
## Correlation tests using Schoenfeld residuals
coxfit=coxph(Surv(remtime,status)~group+logwbc+sex,data=leuk)
sresid=residuals(coxfit,type="schoenfeld")
rankorder=rank(leuk$remtime[leuk$status==1])
cor.test(sresid[,1],rankorder) # group
cor.test(sresid[,2],rankorder) # logwbc
cor.test(sresid[,3],rankorder)$"p.value" # sex
##
## Score tests for time dependent covariates
coxfit=coxph(Surv(remtime,status)~sex,data=leuk)
cox.zph(coxfit)
coxfit=coxph(Surv(remtime,status)~group+logwbc+sex,data=leuk)
testph=cox.zph(coxfit)
testph # chisq values based on score tests
#ggcoxzph(testph) # requires the survminor package
cox.zph(coxfit,log)
tra=function(x) { (x>=7) }
cox.zph(coxfit,tra)
