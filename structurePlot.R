# plot Rscripts for STRUCTURE/Admixture

name_vec<-function(group){
  out=group
  t=as.matrix(table(group))[,1]
  vec<-c()
  vec[1]=ceiling(0.5*t[1])
  for(k in 2:length(t)){
    vec[k]=sum(t[1:k-1])+ceiling(0.5*t[k])
  }
  for(i in 1:length(group)){
    if(i%in%vec==FALSE)
      out[i]=""
  }
  return(out)
}

space_vec=function(group){
  out=c(rep(0,length(group)))
  t=as.matrix(table(group))[,1]
  vec=c()
  for(k in 2:length(t)){
    vec[k-1]=sum(t[1:k-1])
  }
  for(i in vec){
    out[i]=1
  }
  return(out)
}

structurePlot<-function(orderedAdmixtureTable,k,...){
  if(!("Group" %in% names(orderedAdmixtureTable))) 
    stop("Make sure your dataframe contains column Group")
  barplot(t(as.matrix(orderedAdmixtureTable[,1:k])),
          space=space_vec(as.character(orderedAdmixtureTable$Group)),
          border=NA,ylab="Ancestry",
          col=rainbow(k),las=2,
          names.arg = name_vec(as.character(orderedAdmixtureTable$Group)),
          ...
          )
}
