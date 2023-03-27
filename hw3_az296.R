library(survival)
library(rms)
df<-read.csv("C:/Users/andre/Desktop/Cornell/Junior Spring/STSCI4270_Survival_Analysis/Homeworks/Homework 3/gbcs.csv")

#Question 1
## The variables of interest are id,rectime,censrec,size,nodes,grade
df<-df[,c("id","rectime","censrec","size","nodes","grade")]

#Question 1.a
df<-subset(df,grade %in% c(1,3))
df$grade<-as.factor(df$grade)
head(df)

#Question 1.b
m1<-mean(df$rectime[df$grade==1])
m3<-mean(df$rectime[df$grade==3])
## This means that a higher recurrence time is associated with grade 1
## This suggests that grade=1 predicts longer survival than grade=0

#Question 1.c
## KM analysis for grades
kmfit_1.c=npsurv(Surv(rectime,censrec)~grade,data=df) 
## ln[-lnS] plot: grade1 vs grade3
plot(kmfit_1.c,lty=c("solid","dashed"),col=c("black","blue"),
     fun="cloglog",main="ln[-ln S(t)] plot")
legend("topleft",c("grade1","grade3"),
       lty=c("solid","dashed"),col=c("black","blue"))
## Our log-log plots that were plotted from both methods sow that our
## curves do not intersect at any point in time. This suggests that the
## proportional hazards assumption is met for grade

#Question 1.d
## Fit Cox PH on grade only model
coxphfit_1.d=coxph(Surv(rectime,censrec)~grade,data=df)
kmfit_1.d=npsurv(Surv(rectime,censrec)~grade,data=df)
## Observed versus expected plot for grade
survplot(kmfit_1.d,conf.int=FALSE,lty=1,lwd=2,col="black")
newdat=data.frame(grade=c("1","3"))
lines(npsurv(coxphfit_1.d,newdata=newdat),col="red",lty=1,lwd=2)
## Unfortunately, it appears that our residuals do vary across time.
## For grade 1, our residuals are somewhat positive at earlier times
## before becoming somewhat negative in later times. For grade 2, our 
## residuals are somewhat negative before becoming somewhat positive
## at later times. This suggests that grade does not follow the PH assumption

#Question 1.e
## Fit Cox PH size+nodes+grade only model
coxphfit_1.e=coxph(Surv(rectime,censrec)~size+nodes+grade,data=df)
summary(coxphfit_1.e)

#Question 1.f
## Fit Cox PH size+nodes+grade & interactions (grade,size);(grade,nodes) model
coxphfit_1.f=coxph(Surv(rectime,censrec)~size+nodes+grade+grade:size+grade:nodes,data=df)
loglik_1.e=logLik(coxphfit_1.e)
loglik_1.f=logLik(coxphfit_1.f)
-2*(loglik_1.e-loglik_1.f)
## According to our results, our p value is 0.31, this means
## that we fail to reject the null hypothesis. This means that 
## there is no significant difference between our models and that 
## our interaction terms do not add any substantial predictive
## value. 

#Question 1.g
## Score tests for time dependent covariates
coxphfit_1.g=coxph(Surv(rectime,censrec)~size+nodes+grade,data=df)
cox.zph(coxphfit_1.g)
## From our results, it appears that all our covariates
## maintain proportional hazards except for grade. Thus, we 
## can conclude that grade is a time-dependent covariate, also
## meaning that our PH assumption fails for grade.

#Question 1.h
## Fit Cox PH size+nodes (stratified by grade)
coxphfit_1.h=coxph(Surv(rectime,censrec)~size+nodes+strata(grade),data=df)
summary(coxphfit_1.h)
## When our survival model stratifies by grade, we get the coefficient
## 0.071993 for nodes. To relate this to our hazard rate must take the
## exp() of our coefficient.Our exp(coef) ends up being 1.074648
## This means that when we stratify by grade and fit a coxPHmodel on 
## size and nodes, adjusting for size, an increase by 1 unit of nodes is associated with an 
## approximately 1.075 times higher hazard rate. With a p value of 1.46e-08, this effect
## is very significant

#Question 1.i
## Here we specify our slightly modified cplayout function
cplaydf2=function(df1) {
  ## Create data frame, df2, for use with tmerge
  ## in order to generate a data frame in CP layout
  ## Data frame df1 must contain variables: id, rectime, censrec, grade
  ## as well as other variable required for the analysis
  id=df$id
  time=df$rectime
  status=df$censrec
  exposure=df$grade
  etimes=sort(unique(time[status==1])) # ordered unique event times
  n=nrow(df1)
  idc=rep(0,n) # to store number of rows per subject
  tm2=NULL # to store our times for each subject
  for (i in 1:n) { 
    nprior=sum(time[i]>etimes) # number of prior events
    ## We calculate the number of rows needed per id from top to bottom of df$id
    idc[i]=1+nprior # 1+number of prior events
    ## Then we calculate the times to put into the 1+nprior rows below
    if (nprior==0) tm2=c(tm2,time[i])
    if (nprior>0) tm2=c(tm2,etimes[1:nprior],time[i])
  }
  ## We repeat our ids according to nrows specified in idc
  id2=rep(id,idc)
  ## We repeat our exposure according to nrows specified in idc
  ## Our replication per row is exactly the same formula as id replication
  ex2=rep(exposure,idc) 
  df2=data.frame(id=id2,time=tm2,exposure=ex2,exposure.time=as.numeric(ex2)*tm2)
  return(df2)
}
## Now we contruct df2
df2<-cplaydf2(df)
##To construct our (Start,Stop) or Counting Process (CP) format we have more to do
df1<-tmerge(df,df,id,event=event(rectime,censrec))
colnames(df1)[2]="time"
CP_df<-tmerge(df1,df2,id,exposure.tstart=tdc(time,exposure.time,0))
CP_df$exposure.tstop<-as.numeric(CP_df$grade)*CP_df$tstop
##We have finally constructed our CP layout and stored it in CP_df
head(CP_df)

#Question 1.j
## We use our CP Data file to code our Extended model with size,
## nodes, grade, and grade:time as predictors.
coxfit_1.j=coxph(Surv(tstart,tstop,event)~size+nodes+grade+exposure.tstop,data=CP_df)
## Our "exposure" variable is really just grade if you check my table
summary(coxfit_1.j)

#Question 1.k

#Part 1
##See Handwritten Work

#Part 2
##See Handwritten work  
