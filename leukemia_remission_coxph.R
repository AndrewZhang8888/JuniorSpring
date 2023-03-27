## Cox PH analysis of leukemia remission data
## Likelihood-ratio tests
library(survival)
## Read in leukemia remission data
leuk=read.csv("../Data/leukemia_remission.csv")
leuk$group=as.factor(leuk$group)
## Fit Cox PH baseline model
coxphfit0=coxph(Surv(remtime,status)~1,data=leuk)
logLik(coxphfit0)
logLik(coxphfit0)[1]
## Fit Cox PH group only model
coxphfit1=coxph(Surv(remtime,status)~group,data=leuk)
## Plots of survival curves for Treatment and Placebo groups ignoring logWBC
coxphfit.group=with(leuk,data.frame(group=c("0","1"),logwbc=rep(mean(logwbc),2)))
plot(survfit(coxphfit1,newdata=coxphfit.group),
     lwd=2,col=c("black","green"),
     xlab="Remission Time in Weeks",ylab="Survival Probability")
legend("topright",legend=c("Treatment","Placebo"),
       lwd=2,col=c("black","green"))
## Fit Cox PH group+logWBC only model
coxphfit2=coxph(Surv(remtime,status)~group+logwbc,data=leuk)
anova(coxphfit1,coxphfit2)
## Plots of adjusted survival curves for Treatment and Placebo groups 
## controlling for logWBC
plot(survfit(coxphfit2,newdata=coxphfit.group),col=c("black","green"),
     xlab="Remission Time in Weeks",ylab="Survival Probability")
legend("topright",legend=c("Treatment","Placebo"),lwd=2,col=c("black","green"))
##
coxphfit3=coxph(Surv(remtime,status)~group*logwbc,data=leuk)
summary(coxphfit3)
anova(coxphfit2,coxphfit3)
