cplayout=function(df1) {
  ## Create data frame, df2, for use with tmerge
  ## in order to generate a data frame in CP layout
  ## Data frame df1 must contain variables: id, time, status, exposure
  ## as well as other variable required for the analysis
  id=df1$id
  time=df1$time
  status=df1$status
  exposure=df1$exposure
  etimes=sort(unique(time[status==1])) # ordered unique event times
  n=nrow(df1)
  idc=rep(0,n) # to store number of rows per subject
  tm2=NULL # to store times for each subject
  for (i in 1:n) { 
    nprior=sum(time[i]>etimes) # number of prior events
    idc[i]=1+nprior # 1+number of prior events
    if (nprior==0) tm2=c(tm2,time[i])
    if (nprior>0) tm2=c(tm2,etimes[1:nprior],time[i])
  }
  id2=rep(id,idc)
  ex2=rep(exposure,idc)
  df2=data.frame(id=id2,time=tm2,exposure.time=ex2*tm2)
  df2
}
