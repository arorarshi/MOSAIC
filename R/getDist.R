###########################
# getDist
###########################

#computes weights for a given data, standalone
getWeights<-function(xx, out,cv=FALSE, train.snames=NULL){
  #add a check for missing rownames
  if(is.null(rownames(xx))){stop("provide rownames in getDist")}
  #how can we compute faster HRs, apply is faster than for

  if(cv==FALSE){
    #always intersect with out
    out = na.omit(out)
    inter = intersect(names(out),rownames(xx))
    xx = xx[inter,,drop=FALSE]
    out = out[inter]

    #calculate survobj upfront
    ll<-list()
    for(i in 1:ncol(xx)){
      if(length(unique(out))==2){
        ll[[i]]<-summary(glm(out ~ xx[,i], family="binomial"))$coefficients
      }
    }
    if(length(unique(out)) == 2){
      lodds<-unlist(lapply(ll, function(x) x[2,1]))}

    if(length(unique(out)) > 2){

      lodds<-get.probRwts(xx, out)}

    xx.wt <- t(t(xx) * sqrt(abs(lodds)))
    return(xx.wt)
  }

  if(cv==TRUE){
    out=na.omit(out)
    inter = intersect(names(out), intersect(rownames(xx),train.snames))
    xx.train = xx[inter,,drop=FALSE]
    out = out[inter]

    ll<-list()
    for(i in 1:ncol(xx)){
      if(length(unique(out)) == 2){
        ll[[i]]<-summary(glm(out ~ xx.train[,i], family="binomial"))$coefficients
      }
    }

    if(length(unique(out)) == 2){
      lodds<-unlist(lapply(ll, function(x) x[2,1]))}

    if(length(unique(out)) > 2){
      lodds<-get.probRwts(xx.train, out)}

    #on features which are columns, same in all and train
    xx.train.wt <- t(t(xx.train) * sqrt(abs(lodds)))
    xx.wt <- t(t(xx) * sqrt(abs(lodds)))
    wt.mat<-list(all=xx.wt, train=xx.train.wt)
    return(wt.mat)
  }

}

getUnionDist<-function(rnames,dat, type=NULL){
  #compute distance across union of genes/features

  dist.dat<-list()

  for (i in 1:length(dat)){
    m <- dat[[i]][intersect(rownames(dat[[i]]),rnames),,drop=FALSE]
    #if there are missing samples, add row of NAs for it
    if (length(intersect(rownames(m),rnames)) != length(rnames)){
      m.na<-matrix(NA,nrow=length(setdiff(rnames,rownames(m))),ncol=ncol(m))
      rownames(m.na) = setdiff(rnames,rownames(m))
      m<-rbind(m,m.na)
      m<-m[rnames,]
    }
    if(!(is.null(type))){
      if(type=="mut"){dist.dat[[i]] = dist_wtbinary(m)
      rownames(dist.dat[[i]]) = colnames(dist.dat[[i]]) = rownames(m)}}

    if(is.null(type)){
      m2 = m*m
      m2.ss = sum(m2, na.rm=T)
      m.tr = m/sqrt(m2.ss)
      dist.dat[[i]] = as.matrix(dist(m.tr, method="euclidean"))
    }
  }
  return(dist.dat)
}

get.probRwts<-function(xx, out){

  #sometimes unique(out) is not in the same order as table(out)
  uout<-names(table(out))

  bigmu = apply(xx,2,function(x) mean(x, na.rm=T))
  bigsd = apply(xx,2, function(x) sd(x, na.rm=T))

  matc = muc = sigc = loglik = list()
  for(i in 1:length(uout)){
    matc[[i]] = xx[which(out==uout[i]),]
    muc[[i]] = apply(matc[[i]],2, function(x) mean(x, na.rm=T))
    sigc[[i]] = apply(matc[[i]],2,function(x) sd(x, na.rm=T))

  }


  bigl<-list()
  for(i in 1:ncol(xx)){

    fc = fcbig = fl = list()
    for(j in 1:length(uout)){

      fc[[j]] =  sapply(matc[[j]][,i], function(x) dnorm(x, mean=muc[[j]][i], sd = sigc[[j]][i], log=TRUE) )
      fcbig[[j]] =  sapply(matc[[j]][,i], function(x) dnorm(x, mean=bigmu[i], sd = bigsd[i], log=TRUE) )
      #if sd is close to zero, no variance, likelihhood can be Inf.
      #replace Inf values with NA
      #take absolute before taking log. some likelihood can be positive?
      fl[[j]]<-abs(log(abs(sum(fc[[j]], na.rm=T)/sum(fcbig[[j]], na.rm=T))))
      fl[[j]][which(is.infinite(fl[[j]]))] = NA
    }
    bigl[[i]]<-fl
  }

  wt.mat <- matrix(unlist(bigl), ncol = length(uout), byrow = TRUE)
  wts =  apply(wt.mat, 1, function(y) max(y, na.rm=T))

  return(wts)
}


getDist<-function(datasets,out,cv=FALSE,train.snames=NULL,type=NULL){
  # add other checks
  if (is.list(datasets) == F)
    stop("datasets must be a list")

  #if out has no events stop
  if (length(unique(out))==1)
    stop("Need more than one category in outcome")

  #convert everything to numeric
  dat<-lapply(datasets, function(x) as.data.frame(plyr::aaply(x,1,as.numeric,.drop=FALSE)) )
  rnames<-unlist(lapply(datasets, function(x) rownames(x)))
  rnames<-unique(rnames)

  if(is.null(rnames))
    stop("rowanmes=NULL, add sample names to matrix of datasets list object")

  dat.wt<-lapply(dat, function(x) getWeights(x,out,cv,train.snames))

  if(cv==TRUE){
    dat.all.wt<-lapply(dat.wt, function(x) x$all)
    dat.train.wt<-lapply(dat.wt, function(x) x$train)

    dat.all.dist<-getUnionDist(rnames, dat.all.wt,type)
    dat.train.dist<-getUnionDist(train.snames, dat.train.wt, type)
    dat.dist<-list(train=dat.train.dist, all=dat.all.dist)
    return(dat.dist)
  }

  if(cv==FALSE){
    dat.dist<-getUnionDist(rnames, dat.wt, type)
    return(dat.dist)
  }

}

