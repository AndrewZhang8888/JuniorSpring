## VA Lung Cancer Study 
## See e.g. Kalbfleisch and Prentice, Section 4.5
library(survival)
library(rms)
vet=read.csv("C:/Users/andre/Desktop/Cornell/Junior Spring/STSCI4270_Survival_Analysis/Code Data/veteran.csv")
#library(MASS); data(VA)
head(vet)
vet$trt=factor(vet$trt)
vet$celltype=factor(vet$celltype)
levels(vet$celltype)
levels(vet$celltype)=c("squamous","large","adeno","small") # order as in KK
vet$prior=vet$prior/10
coxfit=coxph(Surv(time,status)~trt+celltype+karno+diagtime+age+prior,data=vet)
summary(coxfit)
## Perform score tests for PH assumption
cox.zph(coxfit)
## create 8-level celltype-PSbin factor
vet$PSbin=as.numeric(vet$karno>=60)
vet$z=interaction(vet$celltype,vet$PSbin)
levels(vet$z)
table(vet$z)
## confine attention to treatment and age as in KK
coxfit=coxph(Surv(time,status)~trt+age+strata(z),data=vet)
summary(coxfit)
## fit interaction model and perform LR test
coxfitI=coxph(Surv(time,status)~trt+age+z*trt+z:age+strata(z),data=vet)
anova(coxfit,coxfitI)
summary(coxfitI)
## Extract coefficient estimates and variance matrix
b=as.numeric(coxfitI$coef)
i=(1:length(b))[!is.na(b)] # remove indices of NAs
b=b[i]
V=coxfitI$var[i,i] # variance-covariance matrix
## Survival curves for the x.n stratum adjusted for age
stratum.curves=function(x.n) {
  newdat=data.frame(z=rep(x.n,2),trt=c("1","2"),age=rep(mean(vet$age),2))
  survplot(npsurv(coxfitI,newdata=newdat))
  title(main=paste("Survival curves for",x.n,"celltype",sep=" "),cex.main=1.0)
  mtext(side=3,line=0.25,at=0,adj=-0.5,"(adjusted for age)",cex=1.0)
}
## Wald test for b[1]+b[3]
b13=b[1]+b[3]
v13=V[1,1]+V[3,3]+2*V[1,3]
z=b13/sqrt(v13)
1-pchisq(z^2,1)
## Function to calculate Wald p-value for linear combination of coefficients
wald=function(x) {
  v=t(x)%*%V%*%x
  se=sqrt(v)
  z=t(x)%*%b/se
  1-pchisq(z^2,1)
}
## 
stratum.curves("large.0")
xlarge.0=c(1,0,1,rep(0,(length(i)-3)))
wald(xlarge.0)
vet[vet$z=="large.0",]
##
stratum.curves("adeno.0")
xadeno.0=c(1,0,0,1,rep(0,(length(i)-4)))
wald(xadeno.0)
##

