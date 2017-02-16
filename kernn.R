source("/Users/peterfreeman/PROJECTS/PHOTO_Z/CODE/generic.R")

kernn = function(COVARIATE.SHIFT,z.train.L,z.val.L,z.test.L,
                 dist.val.L_train.L,dist.val.U_train.L,
                 dist.test.L_train.L,dist.test.U_train.L,
                 weights.train.L,weights.val.L,weights.test.L,
                 nn.min=2,nn.max=70,nn.num=10,bw.min=1.e-7,bw.max=5.e-3,
                 bw.num=15,z.min=0,z.max=1,z.num=201,boot.num=400)
{
  print("kernn")

  bw.vec   = seq(log(bw.min),log(bw.max),length.out=bw.num)
  bw.vec   = exp(bw.vec)
  n.neigh  = round(seq(nn.min,nn.max,length.out=nn.num))
  loss     = array(NA,dim=c(length(bw.vec),length(n.neigh)))

  for ( ii in 1:length(bw.vec) ) {
    for ( jj in 1:length(n.neigh) ) {
      cat(".")
      if ( COVARIATE.SHIFT == FALSE ) {
        loss[ii,jj] = error_kernn(n.neigh[jj],z.num,bw.vec[ii],
                                  z.min,z.max,z.train.L,dist.val.L_train.L,
                                  z.val.L,boot=F)$mean
      } else {
        loss[ii,jj] = error_kernn_cs(n.neigh[jj],z.num,bw.vec[ii],
                                     z.min,z.max,z.train.L,weights.train.L,
                                     dist.val.U_train.L,dist.val.L_train.L,
                                     z.val.L,weights.val.L,boot=F)$mean
      }
    }
    cat(" ",round(ii/length(bw.vec),3),"\n")
  }

  point.min = which(loss == min(loss,na.rm=T), arr.ind = TRUE)

  best.bw = (bw.vec)[point.min[1]]
  best.nn = (n.neigh)[point.min[2]]

  pred_kernn = function(dist)
  {
    return(
      t(apply(dist,1,function(xx) {
          near = sort(xx,index.return=T)$ix[1:best.nn]
          if ( COVARIATE.SHIFT == FALSE ) {
            den.obj = cden_kernn(z.train.L[near],z.num,best.bw,z.min,z.max)
          } else {
            den.obj = cden_kernn_cs(z.train.L[near],z.num,best.bw,z.min,z.max,
                                    weights.train.L[near])
          }
          return(den.obj$means)
        }))
    )
  }

  pred.test.L = pred_kernn(dist.test.L_train.L)
  pred.test.U = pred_kernn(dist.test.U_train.L)
  pred.val.L  = pred_kernn(dist.val.L_train.L)
  pred.val.U  = pred_kernn(dist.val.U_train.L)

  final.loss  = error_final_generic(pred.test.L,pred.test.U,
                                    z.test.L,weights.test.L,boot=boot.num)

  return(list(pred.test.L=pred.test.L,pred.test.U=pred.test.U,
              pred.val.L=pred.val.L,pred.val.U=pred.val.U,
              final.loss=final.loss,best.bw=best.bw,best.nn=best.nn))
}

cden_kernn = function(z.train.nearest,n.bins=1000,bw,z.min,z.max)
{
  # calculates estimated density for a single observation
  # weights are previously calculated weights, based on unlabeled data,
  # one weight for each traning sample
  bins.medium = seq(z.min,z.max,length.out=n.bins)
  estimates   = apply(as.matrix(bins.medium),1,function(xx) {
                  weights.final=
                    exp(-abs(xx-z.train.nearest)^2/(4*bw))/sqrt(pi*4*bw)
                  return(sum(weights.final)/length(weights.final))
                })
  output               = NULL
  output$means         = estimates
  output$bin.intervals = bins.medium
  return(output)
}

# This version is without the cVec argument.
cden_kernn_cs = function(z.train.nearest,n.bins=1000,bw,
                         z.min=0,z.max=1,weights.train)
{
  # calculates estimated density for a single observation
  # estimate density of z given z's of the nearest neighbour
  bins.medium = seq(z.min,z.max,length.out=n.bins)
  if ( all(weights.train==0) ) weights.train=rep(1,length(weights.train))
  sum.weights = sum(weights.train)
  estimates   = apply(as.matrix(bins.medium),1,function(xx) {
                  weights.final=
                    weights.train*exp(-(xx-z.train.nearest)^2/(4*bw))/
                      sqrt(pi*4*bw)/sum.weights
                  return(sum(weights.final))
                })
  output                = NULL
  output$means          = estimates
  output$bins.intervals = bins.medium
  return(output)
}

error_kernn = function(n.neigh,n.bins,bw.opt,z.min,
                     z.max,z.train.L,dist.val.L_train.L,
                     z.val.L,boot=F)
{
  z.grid        = seq(z.min,z.max,length.out=n.bins)
  output        = NULL
  output$mean   = NULL
  output$seBoot = NULL

  pred = t( apply(dist.val.L_train.L,1,function(xx) {
              nearest = sort(xx,index.return=T)$ix[1:n.neigh]
              den.obj = cden_kernn(z.train.L[nearest],n.bins,bw.opt,z.min,z.max)
              return(den.obj$means)
          }))
  s.square  = mean(colMeans(pred^2))
 
  n           = length(z.val.L)
  pred.obs    = apply(as.matrix(1:n),1,function(xx) {
                  index=which.min(abs(z.val.L[xx]-z.grid))
                  return(pred[xx,index])
                })
  output$mean = 1/2*s.square-mean(pred.obs)
  if ( boot==F ) return(output);

  # Bootstrap
  meanBoot = apply(as.matrix(1:boot),1,function(xx){
               sample.boot = sample(1:n,replace=T)
               pred.boot   = pred[sample.boot,]
               z.val.boot  = z.val.L[sample.boot]
               s.square    = mean(colMeans(pred.boot^2))
               pred.obs    = apply(as.matrix(1:n),1,function(xx) {
                               index=which.min(abs(z.val.boot[xx]-z.grid))
                               return(pred.boot[xx,index])
                             })
               return(1/2*s.square-mean(pred.obs))
             })
  output$seBoot = sqrt(var(meanBoot))
  return(output)
}

error_kernn_cs = function(n.neigh,n.bins,bw.opt,z.min,z.max,z.train.L,
                          weights.train.L,dist.val.U_train.L,
                          dist.val.L_train.L,z.val.L,weight.z.val.L,boot=F)
{
  z.grid        = seq(z.min,z.max,length.out=n.bins)
  output        = NULL
  output$mean   = NULL
  output$seBoot = NULL

  pred.U = t(apply(dist.val.U_train.L,1,function(xx) {
               nearest = sort(xx,index.return=T)$ix[1:n.neigh]
               den.obj = cden_kernn_cs(z.train.L[nearest],n.bins,bw.opt,z.min,
                                       z.max,weights.train.L[nearest])
               return(den.obj$means)
             }))

  pred.L = t(apply(dist.val.L_train.L,1,function(xx) {
               nearest = sort(xx,index.return=T)$ix[1:n.neigh]
               den.obj = cden_kernn_cs(z.train.L[nearest],n.bins,bw.opt,z.min,
                                       z.max,weights.train.L[nearest])
               return(den.obj$means)
             }))

  s.square  = mean(colMeans(pred.U^2))
  n         = length(z.val.L)
  pred.obs  = apply(as.matrix(1:n),1,function(xx) {
                index = which.min(abs(z.val.L[xx]-z.grid))
                return(pred.L[xx,index])
              })

  output$mean = 1/2*s.square-mean(pred.obs*weight.z.val.L)
  if ( boot==F ) return(output)

  # Bootstrap
  meanBoot = apply(as.matrix(1:boot),1,function(xx) {
               sample.boot          = sample(1:n,replace=T)
               pred.boot.L          = pred.L[sample.boot,]
               pred.boot.U          = pred.U[sample.boot,]
               z.val.boot.L         = z.val.L[sample.boot]
               weights.val.boot.L = weight.z.val.L[sample.boot]
               s.square             = mean(colMeans(pred.boot.U^2))
               pred.obs = apply(as.matrix(1:n),1,function(xx) {
                            index = which.min(abs(z.val.boot.L[xx]-z.grid))
                            return(pred.boot.L[xx,index])
                          })
               return(1/2*s.square-mean(pred.obs*weights.val.boot.L))
             })
  output$seBoot = sqrt(var(meanBoot))
  return(output)
}

