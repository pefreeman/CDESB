
w = which(search()=="package:pracma")
if ( length(w) > 0 ) detach("package:pracma",unload=TRUE)

library(pdist)
library(gplots)
library(fields)

random.labeler = function(n.total,n.L,cov.spec,cov.train.L,
                          cov.train.U,best.nn,min.U)
{
  random.L = sample(1:n.total)
  sample.L = rep(NA,n.L)
  ii = 1
  for ( jj in 1:n.total ) {
    index = random.L[jj]
    dist.L = as.matrix(pdist(cov.spec[index,,drop=F],cov.train.L))
    dist.U = as.matrix(pdist(cov.spec[index,,drop=F],cov.train.U))
    beta.hat = estimate_weights_knn(dist.L,dist.U,best.nn)
    if ( beta.hat*best.nn*(dim(dist.U)[2])/(dim(dist.L)[2]) < min.U ) next;
    sample.L[ii] = index
    if ( ii%%1000 == 0 ) cat("Chose ",ii," out of ",jj,"\n")
    ii = ii+1
    if ( ii > n.L ) break;
  }
  if ( jj == n.total ) stop("There are insufficient labeled data.")
  return(sample.L)
}

compute_cde = function(cov.spec,cov.photo,redshift,n.L=15000,n.U=15000,min.U=1,
                       train.frac=7/15,val.frac=3/15,iseed=401,
                       COVARIATE.SHIFT=FALSE,SERIES=TRUE,KERNN=TRUE,NN=TRUE,
                       USE.Z.BOUNDS=FALSE,Z.MIN.SAMPLE=NULL,Z.MAX.SAMPLE=NULL,
                       NN.MIN=2,NN.MAX=30,NN.NUM=15,SHOW.PLOT=FALSE,
                       OUTPUT.DIR="./",BEST.NN.FIX=-9)
{
  set.seed(iseed)

  covariates = rbind(cov.spec,cov.photo)
  covariates = scale(covariates)  
  cov.spec   = covariates[1:nrow(cov.spec),]
  cov.photo  = covariates[-c(1:nrow(cov.spec)),]
  rm(covariates)

  # Prepare U sample
  sample.U    = sample(1:nrow(cov.photo),n.U,replace=F)
  cov.U       = cov.photo[sample.U,]
  n.train.U   = round(train.frac*length(sample.U))
  n.val.U     = round(val.frac*length(sample.U))
  n.test.U    = n.U - n.train.U - n.val.U
  ran.perm.U  = sample(1:nrow(cov.U))
  cov.train.U = as.matrix(cov.U[ran.perm.U[1:n.train.U],])
  cov.val.U   = as.matrix(cov.U[ran.perm.U[(n.train.U+1):(n.train.U+n.val.U)],])
  cov.test.U  = as.matrix(cov.U[ran.perm.U[(n.train.U+n.val.U+1):n.U],])

  # Prepare L sample
  sample.L    = sample(1:nrow(cov.spec),n.L,replace=F)
  cov.L       = cov.spec[sample.L,]
  n.train.L   = round(train.frac*length(sample.L))
  n.val.L     = round(val.frac*length(sample.L))
  n.test.L    = n.L - n.train.L - n.val.L
  ran.perm.L  = sample(1:nrow(cov.L))
  cov.train.L = as.matrix(cov.L[ran.perm.L[1:n.train.L],])
  cov.val.L   = as.matrix(cov.L[ran.perm.L[(n.train.L+1):(n.train.L+n.val.L)],])
  cov.test.L  = as.matrix(cov.L[ran.perm.L[(n.train.L+n.val.L+1):n.L],])

  dist.train.L_train.L = as.matrix(pdist(cov.train.L,cov.train.L))
  dist.val.L_train.L   = as.matrix(pdist(cov.val.L,cov.train.L))
  dist.val.U_train.L   = as.matrix(pdist(cov.val.U,cov.train.L))
  dist.test.L_train.L  = as.matrix(pdist(cov.test.L,cov.train.L))
  dist.test.U_train.L  = as.matrix(pdist(cov.test.U,cov.train.L))
  
  if ( COVARIATE.SHIFT == TRUE ) {
    if ( BEST.NN.FIX < 0 ) {
      n.neigh = round(seq(NN.MIN,NN.MAX,length.out=NN.NUM))
      loss    = rep(NA,length(n.neigh))
      dist.val.L_train.U = as.matrix(pdist(cov.val.L,cov.train.U))
      dist.val.U_train.U = as.matrix(pdist(cov.val.U,cov.train.U))
      for ( ii in 1:length(n.neigh) ) {
        cat(round(ii/length(n.neigh),2)," ")
        ratio.L.val = estimate_weights_knn(dist.val.L_train.L,
                      dist.val.L_train.U,n.neigh[ii])
        ratio.U.val = estimate_weights_knn(dist.val.U_train.L,
                      dist.val.U_train.U,n.neigh[ii])
        loss[ii]    = estimate_error_weights(ratio.L.val,ratio.U.val,se=F)$mean
      }
      cat("\n")
      rm(dist.val.L_train.U,dist.val.U_train.U)
      if ( SHOW.PLOT == TRUE ) plot(n.neigh,loss,pch=18)
      best.nn  = n.neigh[which.min(loss)]
      cat("The optimal number of nearest neighbors is ",best.nn,"\n")
    } else {
      best.nn  = BEST.NN.FIX
      cat("The number of nearest neighbors is set to be ",best.nn,"\n")
    }
    sample.L    = random.labeler(nrow(cov.spec),n.L,cov.spec,cov.train.L,
                                 cov.train.U,best.nn,min.U)
    cov.L       = cov.spec[sample.L,]
    n.train.L   = round(train.frac*length(sample.L))
    n.val.L     = round(val.frac*length(sample.L))
    n.test.L    = n.L - n.train.L - n.val.L
    ran.perm.L  = sample(1:nrow(cov.L))
    cov.train.L = as.matrix(cov.L[ran.perm.L[1:n.train.L],])
    cov.val.L   = as.matrix(cov.L[ran.perm.L[(n.train.L+1):(n.train.L+n.val.L)],])
    cov.test.L  = as.matrix(cov.L[ran.perm.L[(n.train.L+n.val.L+1):n.L],])
    
    dist.train.L_train.L = as.matrix(pdist(cov.train.L,cov.train.L))
    dist.val.L_train.L   = as.matrix(pdist(cov.val.L,cov.train.L))
    dist.val.U_train.L   = as.matrix(pdist(cov.val.U,cov.train.L))
    dist.test.L_train.L  = as.matrix(pdist(cov.test.L,cov.train.L))
    dist.test.U_train.L  = as.matrix(pdist(cov.test.U,cov.train.L))
  }

  z.L.orig    = redshift[sample.L]
  if ( USE.Z.BOUNDS == TRUE ) {
    z.L       = (z.L.orig-Z.MIN.SAMPLE)/(Z.MAX.SAMPLE-Z.MIN.SAMPLE)
  } else {
    z.L       = (z.L.orig-min(z.L.orig))/(max(z.L.orig)-min(z.L.orig))
  }
  z.train.L   = z.L[ran.perm.L[1:n.train.L]]
  z.val.L     = z.L[ran.perm.L[(n.train.L+1):(n.train.L+n.val.L)]]
  z.test.L    = z.L[ran.perm.L[(n.train.L+n.val.L+1):n.L]]

  z.train.L.orig = z.L.orig[ran.perm.L[1:n.train.L]]
  z.val.L.orig   = z.L.orig[ran.perm.L[(n.train.L+1):(n.train.L+n.val.L)]]
  z.test.L.orig  = z.L.orig[ran.perm.L[(n.train.L+n.val.L+1):n.L]]

  ############################################
  ############ Estimate Weights ############## 
  ############################################

  weights.train.L = rep(1,n.train.L)
  weights.val.L   = rep(1,n.val.L)
  weights.test.L  = rep(1,n.test.L)

  weights.test.U  = rep(1,n.test.U)

  if ( COVARIATE.SHIFT == TRUE ) {
    weights.train.L = estimate_weights_knn(
      as.matrix(pdist(cov.train.L[,,drop=F],cov.train.L[,,drop=F])),
      as.matrix(pdist(cov.train.L[,,drop=F],cov.train.U[,,drop=F])),
      best.nn)
    weights.test.L  = estimate_weights_knn(
      as.matrix(pdist(cov.test.L[,,drop=F],cov.train.L[,,drop=F])),
      as.matrix(pdist(cov.test.L[,,drop=F],cov.train.U[,,drop=F])),
      best.nn)
    weights.val.L   = estimate_weights_knn(
      as.matrix(pdist(cov.val.L[,,drop=F],cov.train.L[,,drop=F])),
      as.matrix(pdist(cov.val.L[,,drop=F],cov.train.U[,,drop=F])),
      best.nn)

    weights.test.U  = estimate_weights_knn(
      as.matrix(pdist(cov.test.U[,,drop=F],cov.train.L[,,drop=F])),
      as.matrix(pdist(cov.test.U[,,drop=F],cov.train.U[,,drop=F])),
      best.nn)
  }

  ###################################

  if ( SERIES == TRUE ) {
    series = series(COVARIATE.SHIFT,z.train.L,z.val.L,z.test.L,
                    dist.train.L_train.L,dist.val.L_train.L,dist.val.U_train.L,
                    dist.test.L_train.L,dist.test.U_train.L,
                    weights.val.L,weights.test.L,
                    show.plot=SHOW.PLOT)
    save(COVARIATE.SHIFT,series,
         z.train.L,z.train.L.orig,z.val.L,z.val.L.orig,z.test.L,z.test.L.orig,
         weights.val.L,weights.test.L,weights.test.U,weights.train.L,
         cov.train.L,sample.L,sample.U,ran.perm.L,ran.perm.U,
         n.train.L,n.val.L,n.test.L,n.train.U,n.val.U,n.test.U,
         file=paste(OUTPUT.DIR,"cde_series.Rdata",sep=""))
  }

  ###################################

  if ( KERNN == TRUE ) {
    kernn = kernn(COVARIATE.SHIFT,z.train.L,z.val.L,z.test.L,
                  dist.val.L_train.L,dist.val.U_train.L,
                  dist.test.L_train.L,dist.test.U_train.L,
                  weights.train.L,weights.val.L,weights.test.L,
                  nn.min=2,nn.max=30,nn.num=15)
    save(COVARIATE.SHIFT,kernn,
         z.train.L,z.train.L.orig,z.val.L,z.val.L.orig,z.test.L,z.test.L.orig,
         weights.val.L,weights.test.L,weights.test.U,weights.train.L,
         cov.train.L,sample.L,sample.U,ran.perm.L,ran.perm.U,
         n.train.L,n.val.L,n.test.L,n.train.U,n.val.U,n.test.U,
         file=paste(OUTPUT.DIR,"cde_kernn.Rdata",sep=""))
  }

  ###################################

  if ( NN == TRUE ) {
    nn = nn(COVARIATE.SHIFT,z.train.L,z.val.L,z.test.L,
            dist.val.L_train.L,dist.val.U_train.L,
            dist.test.L_train.L,dist.test.U_train.L,
            weights.train.L,weights.val.L,weights.test.L,
            nn.min=3,nn.max=90,nn.num=30)
    save(COVARIATE.SHIFT,nn,
         z.train.L,z.train.L.orig,z.val.L,z.val.L.orig,z.test.L,z.test.L.orig,
         weights.val.L,weights.test.L,weights.test.U,weights.train.L,
         cov.train.L,sample.L,sample.U,ran.perm.L,ran.perm.U,
         n.train.L,n.val.L,n.test.L,n.train.U,n.val.U,n.test.U,
         file=paste(OUTPUT.DIR,"cde_nn.Rdata",sep=""))
  }

  ###################################

  return(TRUE)
}


