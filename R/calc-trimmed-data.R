calc.trimmed.data<-function(data,fun,...)
{
 #function returns results specified by 'fun' across a set of "affy-trimmed" array data
 #trimmed refers to removal of top and bottom 2% of data
 #fun could be e.g. mean, median, iqr, other
 
 res=NULL #vector for per-array data
 tmp=NULL

 for(i in (1:ncol(data)))
 {
  tmp<-data[,i][data[,i]>quantile(data[,i],probs=0.02)&data[,i]<quantile(data[,i],probs=0.98)]
  res<-c(res,fun(tmp,...))
 }

res

}
