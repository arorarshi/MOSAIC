
##########################
# MOSAIC
###########################

mosaic<-function(combine.dist,out,k, cmd.k=NULL){
  if(is.null(names(out)))
    stop("rowanmes of out can't be NULL")

  if(is.null(rownames(combine.dist)))
    stop("rowanmes of combine.dist can't be NULL")

  if (!requireNamespace("pdist", quietly = TRUE)) {
    stop("pdist package needed for this function to work. Please install it.",
         call. = FALSE)
  }

  inter <- intersect(names(out), rownames(combine.dist))

  #run cmdscale
  combine.dist <- combine.dist[inter,inter]

  if(is.null(cmd.k)){cmd.k =nrow(combine.dist)-1 }
  if(!(is.null(cmd.k))){cmd.k =as.numeric(cmd.k) }

  cmd.combine.dist<-cmdscale(combine.dist,k=cmd.k)
  my.k = as.numeric(k)
  #run kmeans with 100 starts
  fit<-kmeans(cmd.combine.dist,my.k,nstart=100)

  #return fit
  return(fit)
}

#predict test labels on mosaic fitted
predict.test.label<-function(all.cmd,fit,k){
  all.cmd = as.matrix(all.cmd)
  train.snames = names(fit$cluster)
  test.snames = setdiff(rownames(all.cmd),train.snames)

  #where row - samples, col - genes
  centroid = matrix(NA, nrow = k, ncol = ncol(all.cmd))
  for (kk in 1:k) {
    #meaning k clust has one sample. #WARNING #check
    if(is.vector(all.cmd[names(fit$cluster)[which(fit$cluster==kk)],]) & ncol(all.cmd) > 1){
      message(paste0("k=",k, " training cluster has one sample, prediction might be inaccurate"))
      centroid[kk, ]=all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ]
    }

    if (!(is.null(dim(all.cmd[names(fit$cluster)[fit$cluster==kk],])))){
      if(ncol(all.cmd)> 1){centroid[kk, ]=apply(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ], 2, mean)}
    }

    if(ncol(all.cmd)==1){centroid[kk,] = mean(all.cmd[names(fit$cluster)[which(fit$cluster==kk)], ])}
  }

  dist.whole = apply(centroid,1,function(x) as.matrix(pdist::pdist(x,all.cmd)))

  #assign the cluster membership
  dist.labels = apply(dist.whole,1,which.min)
  names(dist.labels) = rownames(all.cmd)
  test.labels = dist.labels[test.snames]

  #is missing a class label via pdist
  if(length(unique(test.labels)) != k){
    message(paste0("k=", k, " was reduced to ", length(unique(test.labels)), " in test label prediction"))}
  #do we really need to relabel?
  #ul<-unique(test.labels)
  #nl<-1:length(ul)
  #tt<-rep(NA, length(test.labels))
  #names(tt) = names(test.labels)
  #for(i in 1:length(ul)){
  # tt[which(test.labels==ul[i])] = nl[i]
  #}
  #test.labels = tt}

  return(list(test.labels = test.labels))
}
