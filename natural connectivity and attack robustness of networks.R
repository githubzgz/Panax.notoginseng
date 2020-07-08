#caclulation of natural connectivity and three types of 'attack' robustness

#caclulation of natural connectivity of network
natural.connectivity<-function (matrix, normalized = TRUE) 
  {
  #matrix: the adjacent matrix of association network
  #normalized: normalize natural connectivity of network with different size
    eig <- eigen(matrix)
    exp.eig <- exp(eig$values)
    nc <- log(mean(exp.eig))
    if (normalized) {
    n <- length(exp.eig)
    nc <- nc/(n - log(n))
    }
    return(nc)
}

#random 'attack': remove nodes from network randomly (30 repeats)
my.ran.del<-function(matrix){
  result<-rep(1,nrow(matrix))
  nrow<-nrow(matrix)
  for(x in 1:30)
  {
    print(x)
    sample<-vector()
    n<-natural.connectivity(matrix)
    matrix<-matrix
    for(y in 1:(nrow-1))
    {
      sample.1<-sample(setdiff(names(matrix),sample),1)
      matrix[,names(matrix)==sample.1]<-0
      matrix[rownames(matrix)==sample.1,]<-0
      sample<-c(sample,sample.1)
      n1<-natural.connectivity(matrix)
      n<-c(n,n1)
    }
    result<-data.frame(result,n)
  }
  result<-result[,-1]
  result$mean<-apply(result,1,mean)
  result$node<-seq(1:nrow)
  result$node<-(result$node-1)/nrow
  return(result)
}

#remove nodes from network in decreasing sequence of degree centrality
degree.decreasing.robust<-function(matrix,node){
  #node: node table with degree centrality of each node 
  node<-node[order(node$degree,decreasing=T),]
  n<-natural.connectivity(matrix)
  nrow<-nrow(matrix)
  name<-row.names(node)
  for(x in 1:(nrow-1)){
    matrix[,names(matrix)==name[x]]<-0
    matrix[rownames(matrix)==name[x],]<-0
    n1<-natural.connectivity(matrix)
    n<-c(n,n1)
  }
  remove<-(seq(1:nrow(matrix))-1)/nrow(matrix)
  n<-data.frame(n,remove)
  n$n<-as.numeric(n$n)
  return(n)
}

#remove nodes from network in decreasing sequence of betweenness centrality
bet.decreasing.robust<-function(matrix,node){
  #node: node table with betweenness centrality of each node
  node<-node[order(node$bet,decreasing=T),]
  n<-natural.connectivity(matrix)
  nrow<-nrow(matrix)
  name<-row.names(node)
  for(x in 1:(nrow-1)){
    matrix[,names(matrix)==name[x]]<-0
    matrix[rownames(matrix)==name[x],]<-0
    n1<-natural.connectivity(matrix)
    n<-c(n,n1)
  }
  remove<-(seq(1:nrow(matrix))-1)/nrow(matrix)
  n<-data.frame(n,remove)
  n$n<-as.numeric(n$n)
  return(n)
}