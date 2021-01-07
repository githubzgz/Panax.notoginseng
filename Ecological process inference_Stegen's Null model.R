#Stegen's Null model -- Inferring the ecological processes under microbial community assembly 
#1.Beta_NTI
Beta_NTI<-function(phylo,comm,beta.reps=999){
  require(picante)
  #phylo:Phylogenetic tree
  #comm: Community feature table with OTU as rows and site as column
  match.phylo.comm = match.phylo.comm(phylo, t(comm))
  beta.mntd.weighted = as.matrix(comdistnt(match.phylo.comm$comm,cophenetic(match.phylo.comm$phy),abundance.weighted=T))
  
  rand.weighted.bMNTD.comp = array(c(-999),dim=c(nrow(match.phylo.comm$comm),nrow(match.phylo.comm$comm),beta.reps))
  for (rep in 1:beta.reps) {
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(match.phylo.comm$comm,taxaShuffle(cophenetic(match.phylo.comm$phy)),abundance.weighted=T,exclude.conspecifics = F))
    print(c(date(),rep))
  }
  
  weighted.bNTI = matrix(c(NA),nrow=nrow(match.phylo.comm$comm),ncol=nrow(match.phylo.comm$comm))
  for(columns in 1:(nrow(match.phylo.comm$comm)-1)) {
    for(rows in (columns+1):nrow(match.phylo.comm$comm)) {
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals)
      rm("rand.vals")
    }
  }
  
  rownames(weighted.bNTI) = rownames(match.phylo.comm$comm);
  colnames(weighted.bNTI) = rownames(match.phylo.comm$comm);
  return(as.dist(weighted.bNTI))
}

#2.RC_bray
raup_crick= function(mycommunity, reps=999){
  require(ecodist)#Using the function: distance()
  #mycommunity: OTU table with sites as rows and OTUs as columns.
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(mycommunity)
  gamma<-ncol(mycommunity)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(mycommunity), row.names(mycommunity)))
  ##make the mycommunity matrix into a new, pres/abs. matrix:
  ceiling(mycommunity/max(mycommunity))->mycommunity.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(mycommunity.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(mycommunity, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(mycommunity)-1)){
    for(null.two in (null.one+1):nrow(mycommunity)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(mycommunity.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(mycommunity[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(mycommunity[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');			
        ##same for com2:
        com2[sample(1:gamma, sum(mycommunity.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(mycommunity[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(mycommunity[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.mycommunity = rbind(com1,com2); # null.mycommunity;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.mycommunity,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(mycommunity[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results<-as.dist(results)
  return(results)
}

#3.Estimation of ecological process
eco.process<-function(bnti,rcbray){
  n<-nrow(bnti)
  n<-n*(n-1)/2
  bnti<-bnti[upper.tri(bnti)]
  rcbray<-rcbray[upper.tri(rcbray)]
  homo.selection<-sum(bnti<(-2))/n*100
  variable.selection<-sum(bnti>2)/n*100
  rcbray.1<-rcbray[abs(bnti)<2]
  homo.disp<-sum(rcbray.1<(-0.95))/n*100
  disp.limi<-sum(rcbray.1>0.95)/n*100
  undominated<-sum(abs(rcbray.1)<=0.95)/n*100
  cat("homo.selection",homo.selection,
      "variable.selection",variable.selection,
      "homo.disp",homo.disp,
      "disp.limitation",disp.limi,
      "undomiate",undominated)
  process<-c("homo.selection","variable.selection",
             "homo.disp","disp.limitation","undomiate")
  value<-c(homo.selection,variable.selection,
           homo.disp,disp.limi,undominated)
  process.d<-data.frame(process=process,value=value)
  colnames(process.d)<-c("process","value")
  return(process.d)
}