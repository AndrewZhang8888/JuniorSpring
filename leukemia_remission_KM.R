## K-M analysis of leukemia remission data
## log-rank test
library(survival)
## Read in leukemia remission data
leuk=read.csv("../Data/leukemia_remission.csv")
leuk$group=as.factor(leuk$group)
## KM analysis of treatment group
leuk.treatment=leuk[leuk$group==0,]
kmfit.treatment=survfit(Surv(remtime,status)~1,data=leuk.treatment)
summary(kmfit.treatment)
plot(kmfit.treatment,xlab="Survival Time in Weeks",ylab="Survival Probabilities")
## Estimates of mean and median survival
abline(h=0.5,lty=3)
print(kmfit.treatment,rmean=max(leuk.treatment$remtime))
## Cumulative hazard
plot(kmfit.treatment,cumhaz=TRUE,
     xlab="Survival Time in Weeks",ylab="Cumulative Hazard")
## KM analysis of placebo group
leuk.placebo=leuk[leuk$group==1,]
kmfit.placebo=survfit(Surv(remtime,status)~1,data=leuk.placebo)
plot(kmfit.placebo,xlab="Survival Time in Weeks",ylab="Survival Probabilities")
abline(h=0.5,lty=3)
print(kmfit.placebo,rmean=max(leuk.placebo$remtime))
## Joint KM analysis
kmfit=survfit(Surv(remtime,status)~group,data=leuk)
plot(kmfit,lty=c("solid","dashed"),col=c("black","blue"),
     xlab="Survival Time in Weeks",ylab="Survival Probabilities")
legend("topright",c("treatment","placebo"),
       lty=c("solid","dashed"),col=c("black","blue"))
## Log Rank Test
survdiff(Surv(remtime,status)~group,data=leuk)
