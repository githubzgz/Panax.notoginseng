#Calculation of niche breadth (B)
##comun: A community table with samples as rows and taxa as columns. 
Com.niche <- function(comun, stats=TRUE){
  require(spaa)
  comun<-comun[,colSums(comun)>0]
  B<-niche.width(comun,method="levins")
  B_com<-1:nrow(comun)
  for(i in 1:nrow(comun)){
    a<-comun[i,]
    a<-a[a>0]
    B_com[i]<-mean(as.numeric(B[,names(a)]))
  }
  return(B_com)
}


#Calculation of relative dispersal potential for bacteria and fungi, respectively
b_dispersal<-data.frame(matrix(NA,nrow=78,ncol=78))
for(i in 1:78){
  colsum<-colSums(b[,1:78])[1]
  for(j in i:78){
    compare<-data.frame(i=b[,i],j=b[,j])
    compare$min<-apply(compare,1,min)
    share<-sum(compare$min)/colsum*100
    b_dispersal[i,j]<-share
    b_dispersal[j,i]<-share
  }
}

f_dispersal<-data.frame(matrix(NA,nrow=78,ncol=78))
for(i in 1:78){
  colsum<-colSums(f[,1:78])[1]
  for(j in i:78){
    compare<-data.frame(i=f[,i],j=f[,j])
    compare$min<-apply(compare,1,min)
    share<-sum(compare$min)/colsum*100
    f_dispersal[i,j]<-share
    f_dispersal[j,i]<-share
  }
}