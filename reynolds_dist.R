#set.seed(44)

#op1<-runif(12)
#pop2<-runif(12)
#df1<-rbind(pop1,pop2)
#reynolds_dist<-reynolds_distance(df1)

reynolds_distance=function(Data) 
{
  popnames=rownames(Data)
  npop=nrow(Data) 
  nloc=ncol(Data) 
  dist=matrix(0,nrow=npop,ncol=npop) 
  for (i in 1:(npop-1)) 
  { 
    for (j in (i+1):npop) 
    { 
      pi=Data[i,] 
      pj=Data[j,] 
      
      dist[i,j]=sqrt(sum( (pi-pj)**2+(pj-pi)**2 )/(2*sum(1-(pi*pj+(1-pi)*(1-pj)))) )
    } 
  } 
  dist=dist+t(dist) 
  rownames(dist)=colnames(dist)=popnames
  return(dist) 
}