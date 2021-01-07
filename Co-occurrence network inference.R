#######################
#Network Inference
#######################
require(CoNetinR)#getPval()
require(RMThreshold)#

#1.The calculation of P value for Spearman correlation and Bray-Curtis distance
#(1)Spearman P value
spearman.p<-matrix(nrow=nrow(all),ncol=nrow(all))

for(i in 1:nrow(all)){
  for(j in 1:i){
    if(i != j){
      spearman.p[i, j] = getPval(as.matrix(all),i,j,method = "spearman",permutandboot = T,renorm = T)
      spearman.p[j, i] = spearman.p[i, j]
    }
    print(c(i,j))
  }
}

#(2)Bray P value
bray.p<-matrix(nrow=nrow(all),ncol=nrow(all))

for(i in 1263:nrow(all)){
  for(j in 1:i){
    if(i != j){
      bray.p[i, j] = getPval(as.matrix(all),i,j,method = "bray",permutandboot = T)
      bray.p[j, i] = bray.p[i, j]
    }
    print(c(i,j))
  }
}

#2.Combination of P values from two measures and correction
q=-2*log(spearman.p)-2*log(bray.p)

merge.p<-matrix(nrow=dim(bray.p)[1],ncol=dim(bray.p)[1])

for(j in 1:dim(bray.p)[1]){
  for(k in 1:dim(bray.p)[1]){
    merge.p[j,k] <- pchisq(q[j,k],4,lower.tail = FALSE)
    #Freedom = 2K (K is the number of measures needed to be merged)
  }
  print(c(j,k))
}

merge.fdr<-matrix(p.adjust(merge.p,method="BH"),nrow=dim(merge.p))
diag(merge.fdr)<-1


#3.The estimation of threshold for Spearman correlations based on the Random Matrix Theory
#(1)Calculating the Spearman correlations
spearman.r<-corr.test(t(all_com),method="spearman",adjust="none",ci=F)#all_com is the combined community feature table including bacterial and fungal OTUs
spearman.cor.r<-spearman.r$r
#(2)The assessment of thresholds
spearman.thresh<-rm.get.threshold(spearman.cor.r,nr.thresholds=20,
                                  wait.seconds = 1,unfold.method="spline")#0.74

