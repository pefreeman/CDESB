source("/Users/peterfreeman/PROJECTS/PHOTO_Z/CODE/generic.R")

nn = function(COVARIATE.SHIFT,z.train.L,z.val.L,z.test.L,
              dist.val.L_train.L,dist.val.U_train.L,
              dist.test.L_train.L,dist.test.U_train.L,
              weights.train.L,weights.val.L,weights.test.L,
              bin.min=5,bin.max=50,bin.num=10,nn.min=5,nn.max=300,nn.num=10,
              z.min=0,z.max=1,z.num=201)
{
  print("nn")

  n.bins  = round(seq(bin.min,bin.max,length.out=bin.num))
  n.neigh = round(seq(nn.min,nn.max,length.out=nn.num))
  loss = array(NA,dim=c(length(n.bins),length(n.neigh)))
  for ( ii in 1:length(n.bins) ) {
    for ( jj in 1:length(n.neigh) ) {
      cat(".")
      if ( COVARIATE.SHIFT == FALSE ) {
        loss[ii,jj] = error_nn(n.neigh[jj],300,n.bins[ii],z.min,z.max,z.train.L,
                               dist.val.L_train.L,z.val.L,boot=F)$mean
      } else {
        loss[ii,jj] = error_nn_cs(n.neigh[jj],300,n.bins[ii],z.min,z.max,z.train.L,
                                  weights.train.L,dist.val.U_train.L,
                                  dist.val.L_train.L,z.val.L,weights.val.L,
                                  boot=F)$mean
      }
    }
    cat(" ",round(ii/length(n.bins),3),"\n")
  }

  point.min = which(loss == min(loss,na.rm=T), arr.ind = TRUE)

  best.bin = (n.bins)[point.min[1]]
  best.nn  = (n.neigh)[point.min[2]]

  z.grid         = seq(z.min,z.max,length.out=z.num)
  bins.intervals = seq(0,1,length.out=best.bin)
  which.closest  = apply(as.matrix(z.grid),1,function(yy) { 
                     which.min(abs(bins.intervals-yy)) 
                   })

  pred_nn = function(dist)
  {
    return(
      t(apply(dist,1,function(xx) {
          near = sort(xx,index.return=T)$ix[1:best.nn]
          if ( COVARIATE.SHIFT == FALSE ) {
            den.obj = cden_nn(z.train.L[near],best.bin,z.min,1)
            means   = den.obj$means[which.closest]
          } else {
            den.obj = cden_nn_cs(z.train.L[near],z.num,z.min,z.max,
                                 weights.train.L[near])
#### Should there be a which.closest here?
            means   = den.obj$means
          }
          return(means)
        }))
    )
  }

  pred.test.L = pred_nn(dist.test.L_train.L)
  pred.test.U = pred_nn(dist.test.U_train.L)
  pred.val.L  = pred_nn(dist.val.L_train.L)
  pred.val.U  = pred_nn(dist.val.U_train.L)

  if ( COVARIATE.SHIFT == FALSE ) {
    final.loss = error_final_generic(pred.test.L,pred.test.U,
                                     z.test.L,weights.test.L,boot=400)
  } else {
    # Why is the following *not* error_final_generic?
    final.loss = error_nn_cs(best.nn,1000,best.bin,z.min,z.max,z.train.L,
                             weights.train.L,dist.test.U_train.L,
                             dist.test.L_train.L,z.test.L,
                             weights.test.L,boot=400)
  }

  return(list(pred.test.L=pred.test.L,pred.test.U=pred.test.U,
              pred.val.L=pred.val.L,pred.val.U=pred.val.U,
              final.loss=final.loss,best.bin=best.bin,best.nn=best.nn,
              which.closest=which.closest))
}

cden_nn = function(z.train.nearest,n.bins,z.min,z.max)
{
  # calculates estimated density for a single observation
  # estimate density of z given z's of the  nearest neighbour
  bin.intervals  = seq(z.min,z.max,length.out=n.bins)
  bin.size       = bin.intervals[2]-bin.intervals[1]
  means          = apply(as.matrix(bin.intervals),1,function(xx) {
                     lower = xx-bin.size/2
                     upper = xx+bin.size/2
                     mean(z.train.nearest>=lower & z.train.nearest<upper)/
                       bin.size
                   })
  output               = NULL
  output$means         = means
  output$bin.intervals = bin.intervals
  return(output)
}

# Version without cVec.
cden_nn_cs = function(z.train.nearest,n.bins,z.min,z.max,weights.train)
{
  # calculates estimated density for a single observation
  # estimate density of z given z's of the  nearest neighbour
  bins.intervals = seq(z.min,z.max,length.out=n.bins)
  bin.size       = bins.intervals[2]-bins.intervals[1]
  estimates      = apply(as.matrix(bins.intervals),1,function(xx) {
                     lower = xx-bin.size/2
                     upper = xx+bin.size/2
                     if ( all(weights.train==0) ) {
                       weights.train = rep(1,length(weights.train))
                     }
                     sum(weights.train[z.train.nearest>=lower &
                       z.train.nearest<upper])/(bin.size*sum(weights.train))
                   })
  output                = NULL
  output$means          = estimates
  output$bins.intervals = bins.intervals
  return(output)
}

error_nn = function(n.neighbors,n.bins,n.bins.opt,z.min,z.max,z.train,
                    dist.val.train,z.val,boot=F)
{
  z.grid        = seq(z.min,z.max,length.out=n.bins)
  bin.intervals = seq(z.min,z.max,length.out=n.bins.opt)
  which.closest = apply(as.matrix(z.grid),1,function(yy)  {
                    which.min(abs(bin.intervals-yy))
                  })

  pred          = t(apply(dist.val.train,1,function(xx) {
                      nearest = sort(xx,index.return=T)$ix[1:n.neighbors]
                      den.obj = cden_nn(z.train[nearest],n.bins.opt,z.min,z.max)
                      means   = den.obj$means[which.closest]
                      return(means)
                  }))

  col.means.comp = colMeans(pred^2)
  s.square       = mean(col.means.comp)
  n              = length(z.val)
  pred.obs       = apply(as.matrix(1:n),1,function(xx) {
                     index=which.min(abs(z.val[xx]-z.grid))
                     return(pred[xx,index])
                   })
  output         = NULL
  output$mean    = 1/2*s.square-mean(pred.obs)
  if ( boot==F ) return(output)

  # Bootstrap
  meanBoot = apply(as.matrix(1:boot),1,function(xx) {
               sample.boot=sample(1:n,replace=T)
               pred.boot=pred[sample.boot,]
               z.val.boot=z.val[sample.boot]
               col.means.comp =colMeans(pred^2)
               s.square=mean(col.means.comp )
               pred.obs=apply(as.matrix(1:n),1,function(xx) {
                          index=which.min(abs(z.val.boot[xx]-z.grid))
                          return(pred.boot[xx,index])
                        })
               return(1/2*s.square-mean(pred.obs))
             })
  output$seBoot = sqrt(var(meanBoot))
  return(output)
}

# Version without cVec.
error_nn_cs = function(n.neigh,n.bins,n.bins.opt,z.min,z.max,z.train.L,
                       weights.train.L,dist.test.U_train.L,
                       dist.test.L_train.L,z.test.L,weights.test.L,boot=F)
{
  z.grid         = seq(z.min,z.max,length.out=n.bins)
  bins.intervals = seq(z.min,z.max,length.out=n.bins.opt)
  which.closest  = apply(as.matrix(z.grid),1,function(yy) {
                     which.min(abs(bins.intervals-yy))
                   })
  output         = NULL
  output$mean    = NULL
  output$seBoot  = NULL

  pred.U = t(apply(dist.test.U_train.L,1,function(xx) {
               nearest = sort(xx,index.return=T)$ix[1:n.neigh]
               den.obj = cden_nn_cs(z.train.L[nearest],n.bins.opt,
                                    z.min,z.max,weights.train.L[nearest])
               means   = den.obj$means[which.closest]
               return(means)
             }))

  pred.L = t(apply(dist.test.L_train.L,1,function(xx) {
               nearest = sort(xx,index.return=T)$ix[1:n.neigh]
               den.obj = cden_nn_cs(z.train.L[nearest],n.bins.opt,
                                    z.min,z.max,weights.train.L[nearest])
               means   = den.obj$means[which.closest]
               return(means)
             }))

  col.means = colMeans(pred.U^2)
  s.square  = mean(col.means)
  n         = length(z.test.L)
  pred.obs  = apply(as.matrix(1:n),1,function(xx) {
                index = which.min(abs(z.test.L[xx]-z.grid))
                return(pred.L[xx,index])
              })

  output$mean = 1/2*s.square-mean(pred.obs*weights.test.L)
  if ( boot==F ) return(output)

  # Bootstrap
  meanBoot = apply(as.matrix(1:boot),1,function(xx) {
               sample.boot           = sample(1:n,replace=T)
               pred.boot.L           = pred.L[sample.boot,]
               pred.boot.U           = pred.U[sample.boot,]
               z.test.boot.L         = z.test.L[sample.boot]
               weights.test.boot.L   = weights.test.L[sample.boot]
               col.means             = colMeans(pred.boot.U^2)
               s.square              = mean(col.means)
               pred.obs = apply(as.matrix(1:n),1,function(xx) {
                            index=which.min(abs(z.test.boot.L[xx]-z.grid))
                            return(pred.boot.L[xx,index])
                          })
               return(1/2*s.square-mean(pred.obs*weights.test.boot.L))
             })
  output$seBoot = sqrt(var(meanBoot))
  return(output)
}

