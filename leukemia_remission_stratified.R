## Leukemia remission data
## Stratified Cox Model
library(survival)
library(rms)
#library(simPH) # required for ggfitStrata
#library(survminer) # required for ggcoxzph()
## Read in leukemia remission data
leuk=read.csv("../Data/leukemia_remission_with_sex.csv")
leuk$group=as.factor(leuk$group)
leuk$sex=as.factor(leuk$sex)
## Full model stratified by Sex
coxfit1=coxph(Surv(remtime,status)~group+logwbc+strata(sex),data=leuk)
loglik1=logLik(coxfit1)
## Reduced model (without group) stratified by Sex 
coxfit0=coxph(Surv(remtime,status)~logwbc+strata(sex),data=leuk)
loglik0=logLik(coxfit0)
## Likelihood-ratio test for group controlling for logWBC and Sex
-2*(loglik0-loglik1)
anova(coxfit0,coxfit1)
## Separate models for males and females
coxfitM=coxph(Surv(remtime,status)~group+logwbc,data=leuk,subset=(sex==0))
loglikM=logLik(coxfitM)
coxfitF=coxph(Surv(remtime,status)~group+logwbc,data=leuk,subset=(sex==1))
loglikF=logLik(coxfitF)
loglikMF=loglikM+loglikF
## Alternative specification of interaction model
coxfitI=coxph(Surv(remtime,status)~group+logwbc+sex:group+sex:logwbc+strata(sex),data=leuk)
loglikI=logLik(coxfitI)
anova(coxfit1,coxfitI)
## Plot of fitted survival curves by sex and treatment adjusted for logWBC
newdat=data.frame(group=rep(c("0","1"),2),logwbc=rep(mean(leuk$logwbc),4),sex=rep(c("0","1"),each=2))
#newdat=with(leuk,data.frame(group=rep(c("0","1"),2),logwbc=rep(mean(leuk$logwbc),4),sex=rep(c("0","1"),each=2)))
survplot(npsurv(coxfit1,newdata=newdat))
##ggfitStrata(npsurv(coxfit1,newdata=newdat))
