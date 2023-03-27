library(survival)
library(survminer)
library(broom)

stime=c(10,15,23,23,30,31,35,52,100)
status=c(1,1,1,1,1,0,1,1,0)

#Question 1.A
myKM=function(stime,status){
  surv<-numeric()
  for (i in 1:length(stime)){
    ni<-length(stime)-i+1
    mi<-status[i]
    si<-(ni-mi)/ni
    if(i==1){
      surv<-append(surv,si)
    }
    else{
      si<-(ni-mi)/ni*surv[i-1]
      surv<-append(surv,si)
    }
  }
  return(surv)
}

#Question 1.B
codeKM<-survfit(Surv(stime,status)~1)
summary(codeKM)
print(myKM(stime,status))

#Our Survival Function Estimates at each time point are the same
#when comparing my manual function and the R Package

#Question 1.C
plot(codeKM, xlab='Time',ylab='Survival Probability')

#Our median survival time is 30 days according to the 0.5 quantile of Survival Probability
#Our 95% CI for median survival time can be found by drawing a horizontal line at the 0.5 quantile:
#[20 days,NA]
#Our upper bound is undefined since the upper bound of survival time stretches off to infinity

#Question 1.D
sum_vec<-numeric()
KM_table<-tidy(codeKM)
for (i in 1:7){
  risk_set=KM_table[i,"n.risk"]
  events=KM_table[i,"n.event"]
  local_sum=events/(risk_set*(risk_set-events))
  sum_vec<-append(sum_vec,local_sum)
}
total_sum=do.call(sum,sum_vec)
var<-0.148^2*total_sum
stderr<-sqrt(var)
z<-1.96
expected<-KM_table[7,"estimate"]
low_bound<-expected-z*stderr
high_bound<-expected+z*stderr

#Our Confidence interval for the probability of survival past one year is 
#[0,0.408529]
#Disclaimer:
#The generated lower bound was given as -0.1122
#However, I've adjusted the lower bound to 0 as a negative 
#survival probability cannot logically be interpreted

#Question 1.E
for (i in 1:7){
  risk_set=KM_table[i,"n.risk"]
  events=KM_table[i,"n.event"]
  local_sum=events/(risk_set*(risk_set-events))
  sum_vec<-append(sum_vec,local_sum)
}
var<-0.805555
log_mean<-log(KM_table[7,"estimate"])
z<-1.96
# Our 95%CI for log(s(t)) is:
# [-3.6687,-0.15039]
# If we exp() both bounds, we get the approximate 95%CI for s(t)
# [0.02551,0.860372]

#Question 1.F

#I calculated the work by hand in an attachment
#I find it much easier to properly do the delta method on paper
#Our final 95%CI interval is:
#[0.00825,0.467646]

#Question 2.A
agvhd<-read.csv("C:/Users/andre/Desktop/Cornell/Junior Spring/STSCI4270_Survival_Analysis/Homeworks/Homework 2/AGVHD.csv")
agvhd$Group<-as.factor(agvhd$Group)
agvhd_surv<-survfit(Surv(Time,Status)~Group,data=agvhd)
ggsurvplot(
  fit=agvhd_surv,
  data=agvhd,
  xlab='Time (days)',
  ylab='Survival Probability',
  xlim=c(0,80),
  legend='bottom',
  legend.title="",
  legend.labs=c("Group0","Group1"),
)

#Question 2.B
survdiff(Surv(Time,Status)~Group,data=agvhd)
#According to our log-rank test, we have a pvalue of 0.02
#This means that our 2 Kaplan-Meier Curves are significantly different

#Question 2.C
agvhd<-agvhd[order(agvhd$Time),]
agvhd$index<-1:nrow(agvhd)
agvhd$n_0<-32
agvhd$n_1<-32
agvhd$m_0<-0
agvhd$m_1<-0
for (i in 1:63){
  if(agvhd[i,c('Group')]==0){
    agvhd[i+1,c('n_0')]=agvhd[i,c('n_0')]-1
    agvhd[i+1,c('n_1')]=agvhd[i,c('n_1')]
    if(agvhd[i,c('Status')]==1){
      agvhd[i,c('m_0')]=1
    }
  }
  else{
    agvhd[i+1,c('n_1')]=agvhd[i,c('n_1')]-1
    agvhd[i+1,c('n_0')]=agvhd[i,c('n_0')]
    if(agvhd[i,c('Status')]==1){
      agvhd[i,c('m_1')]=1
    }
  }
}
agvhd$e_0<-0
agvhd$e_1<-0
for (i in 1:64){
  p_0=agvhd[i,c('n_0')]/(agvhd[i,c('n_1')]+agvhd[i,c('n_0')])
  p_1=agvhd[i,c('n_1')]/(agvhd[i,c('n_1')]+agvhd[i,c('n_0')])
  agvhd[i,c('e_0')]=(agvhd[i,c('m_0')]+agvhd[i,c('m_1')])*p_0
  agvhd[i,c('e_1')]=(agvhd[i,c('m_0')]+agvhd[i,c('m_1')])*p_1
}
O_0<-do.call(sum,agvhd['m_0'])
O_1<-do.call(sum,agvhd['m_1'])
E_0<-do.call(sum,agvhd['e_0'])
E_1<-do.call(sum,agvhd['e_1'])
g0_stat<-(O_0-E_0)^2/E_0
g1_stat<-(O_1-E_1)^2/E_1
log_stat<-g0_stat+g1_stat
p_val<-pchisq(log_stat,df=1,lower.tail=FALSE)

#Our Observed Counts in group 0 is: 5
#Our Observed Counts in group 1 is: 15
#Our Expected Counts in group 0 is: 10.2195
#Our Expected Counts in group 1 is: 9.7805
#Our calculated p-value is 0.01955
#This means that we reject the null hypothesis that the two survival curves are the same

#Question 2.D
cox<-coxph(Surv(Time,Status)~Group,data=agvhd)
summary(cox)

#Our exp(coef) for the group variable is 3.1645
#This means that a person in group_1 is 3.16 times more 
#likely to experience an event or failure at any given
#instant in time. 

#Our exp(-coef) for the group variable is 0.316. This
#means that a person in group0 is 0.316 times as likely
#to experience an event or failure compared to a 
#person in group1 at any given time. 

#Question 2.E
cox1<-coxph(Surv(Time,Status)~Group+Age,data=agvhd)
summary(cox1)

#Our exp(coef) for the group variable is 3.957
#This means that a person in group_1 is 3.957 times as 
#likely to experience an event or failure at any given
#instant in time, when adjusted for age value!

#Our exp(-coef) for the group variable is 0.2527. This
#means that a person in group0 is 0.2527 times as likely
#to experience an event or failure compared to a 
#person in group1 at any given time, when adjusted
#for age value!

#Question 2.F
agvhd$AgeInd<-ifelse(agvhd$Age>=16,1,0)
agvhd$AgeInd<-as.factor(agvhd$AgeInd)
ageind2<-unlist(agvhd['AgeInd'])
cox2<-coxph(Surv(Time,Status)~Group+AgeInd,data=agvhd)
summary(cox2)

#Our exp(coef) for the ageind variable is 6.206
#This means that a person over the age of 16 is 6.206 times as 
#likely to experience an event or failure at any given
#instant in time compared to a person under the age of 16, when adjusted for group value!

#Our exp(-coef) for the ageind variable is 0.1611
#This means that a person under the age of 16 is 0.1611 times as
#likely to experience an event or failure at any given
#instant in time compared to a person of age 16 or older, when adjusted for group value!

#Question 2.G
cox3<-coxph(Surv(Time,Status)~AgeInd+Group,data=agvhd,subset=(AgeInd==1))
ggsurvplot(
  fit=survfit(cox3),
  data=agvhd,
  xlab='Time (days)',
  ylab='Survival Probability',
  xlim=c(0,80),
  legend='bottom',
  legend.title="",
  legend.labs=c("Groups"),
)

#Question 2.H
cox4<-coxph(Surv(Time,Status)~AgeInd+Group,data=agvhd,subset=(Group==1))
ggsurvplot(
  fit=survfit(cox4),
  data=agvhd,
  xlab='Time (days)',
  ylab='Survival Probability',
  xlim=c(0,80),
  legend='bottom',
  legend.title="",
  legend.labs=c("AgeInd"),
)

#Question 2.I
cox5<-coxph(Surv(Time,Status)~AgeInd+Group+Group:AgeInd,data=agvhd,subset=(Group==1))
ggsurvplot(
  fit=survfit(cox5),
  data=agvhd,
  xlab='Time (days)',
  ylab='Survival Probability',
  xlim=c(0,80),
  legend='bottom',
  legend.title="",
  legend.labs=c("AgeInd"),
)

anova(cox5,cox4)
#According to my ANOVA test, the 2 survival curves are identical.
#I guess that by taking a subset of "agvhd$Group==1" we have made 
#aghvd$Group a constant regressor. It is intuitive that 
#a random variable has zero interaction with a constant number,
#meaning that Group:AgeInd is always an interaction effect of 0!
