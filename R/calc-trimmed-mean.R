calc.trimmed.mean <- function(data){
 #returns mean of data, excluding values that are >=98-percentile or <=2-percentile
 #used in affy scaling-normalisatino scheme
 
  res=NULL
  
  res<-mean(data[data>quantile(data,probs=0.02)&data<quantile(data,probs=0.98)])
  
  res
}
 
