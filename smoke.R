library(survival)
id1=100+1:4
tm1=c(2,3,5,8)
st1=c(1,1,0,1)
ex1=c(1,0,0,1)
df1=data.frame(id=id1,time=tm1,status=st1,exposure=ex1)
##
##########################
## Compute CP data layout 
etimes=sort(unique(df1$time[df1$status==1])) # ordered unique event times
n=nrow(df1)
idc=rep(0,n) # to store number of rows per subject
tm2=NULL # to store times for each subject
for (i in 1:n) { 
  nprior=sum(df1$time[i]>etimes) # number of prior events
  idc[i]=1+nprior # 1+number of prior events
  if (nprior==0) tm2=c(tm2,df1$time[i])
  if (nprior>0) tm2=c(tm2,etimes[1:nprior],df1$time[i])
}
id2=rep(df1$id,idc)
ex2=rep(df1$exposure,idc)
df2=data.frame(id=id2,time=tm2,exposure.time=ex2*tm2)
## Use tmerge to compute complete CP data layout
df1=tmerge(df1,df1,id,event=event(time,status))
newdf=tmerge(df1,df2,id,exposure.tstart=tdc(time,exposure.time))
newdf$exposure.tstop=newdf$exposure*newdf$tstop
## End of code for time-dependent covariates
##########################
