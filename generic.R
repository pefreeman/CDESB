
error_final_generic = function(pred.L.test,pred.U.test,z.test.L,
                               weights.test.L,boot=F)
{
  z.grid         = seq(0,1,length.out=dim(pred.L.test)[2])
  col.means.comp = colMeans(pred.U.test^2)
  s.square       = mean(col.means.comp)
  
  n.U            = dim(pred.U.test)[1]
  n.L            = length(z.test.L)
  pred.obs       = apply(as.matrix(1:n.L),1,function(xx) { 
                     index=which.min(abs(z.test.L[xx]-z.grid))
                     return(pred.L.test[xx,index])
                   })
  output         = NULL
  output$mean    = 1/2*s.square-mean(pred.obs*weights.test.L)
  if ( boot==F ) return(output)
  
  # Bootstrap
  meanBoot = apply(as.matrix(1:boot),1,function(xx) {
               sample.boot.L=sample(1:n.L,replace=T)
               sample.boot.U=sample(1:n.U,replace=T)
               pred.boot.L=pred.L.test[sample.boot.L,]
               pred.boot.U=pred.U.test[sample.boot.U,]
               z.test.L.boot=z.test.L[sample.boot.L]
               weights.test.boot.L=weights.test.L[sample.boot.L]
               col.means.comp =colMeans(pred.boot.U^2)
               s.square=mean(col.means.comp)
               pred.obs = apply(as.matrix(1:n.L),1,function(xx) { 
                            index=which.min(abs(z.test.L.boot[xx]-z.grid))
                            return(pred.boot.L[xx,index])
                          })
               return(1/2*s.square-mean(pred.obs*weights.test.boot.L))    
             })
  output$seBoot = sqrt(var(meanBoot))
  return(output)
}

normalize_density = function(bin.size,est,delta=0)
{
  est                      = matrix(est,1,length(est))
  if(all(est<=0)) est      = matrix(1,1,length(est))
  est.thresh               = est
  est.thresh[est.thresh<0] = 0
  
  if ( sum(bin.size*est.thresh) > 1 ) {
    max.dens = max(est)  
    min.dens = 0
    newXi    = (max.dens+min.dens)/2
    eps      = 1
    ii       = 1
    while ( ii <= 1000 ) {
      est.new = apply(as.matrix(est),2,function(xx)max(0,xx-newXi))
      area          = sum(bin.size*est.new)
      eps           = abs(1-area)
      if( eps < 0.0000001) break; # level found
      if( 1 > area) max.dens = newXi
      if( 1 < area) min.dens = newXi
      newXi = (max.dens+min.dens)/2
      ii    = ii+1
    }
    est.new = apply(as.matrix(est),2,function(xx)max(0,xx-newXi))
    
    runs   = rle(est.new>0)
    n.runs = length(runs$values)
    jj     = 1
    area   = lower = upper = NULL
    if ( n.runs > 2 ) {
      for ( ii in 1:n.runs ) {
        if ( runs$values[ii] == FALSE ) next;
        whichMin = 1
        if ( ii > 1 ) whichMin = sum(runs$lengths[1:(ii-1)])
        whichMax = whichMin+runs$lengths[ii]
        lower[jj] = whichMin # lower interval of component
        upper[jj] = whichMax # upper interval of component
        # total area of component
        area[jj]  = sum(bin.size*est.new[whichMin:whichMax]) 
        jj        = jj+1
      }
      
      delta = min(delta,max(area))
      for ( ii in 1:length(area) ) {
        if ( area[ii] < delta ) est.new[lower[ii]:upper[ii]] = 0
      }
      est.new = est.new/(bin.size*sum(est.new))
    }
    return(est.new)
  }
  
  est.new = as.vector(1/bin.size*est.thresh/sum(est.thresh))
  
  runs   = rle(est.new>0)
  n.runs = length(runs$values)
  jj     = 1
  area   = lower = upper = NULL
  if ( n.runs > 2 ) {
    for ( ii in 1:n.runs ) {
      if ( runs$values[ii] == FALSE ) next;
      whichMin = 1
      if ( ii > 1 ) whichMin = sum(runs$lengths[1:(ii-1)])
      whichMax  = whichMin+runs$lengths[ii]
      lower[jj] = whichMin # lower interval of component
      upper[jj] = whichMax # upper interval of component
      # total area of component
      area[jj]  = sum(bin.size*est.new[whichMin:whichMax]) 
      jj        = jj+1
    }
    delta = min(delta,max(area))
    for ( ii in 1:length(area) ) {
      if ( area[ii] < delta ) est.new[lower[ii]:upper[ii]] = 0
    }
    est.new = est.new/(bin.size*sum(est.new))
  }
  return(est.new)
}

estimate_weights_knn = function(dist.X.test.train.L,dist.X.test.train.U,n.neigh)
{
  # outputs the estimate at the points XTest
  n.L    = dim(dist.X.test.train.L)[2]
  n.U    = dim(dist.X.test.train.U)[2]
  n.test = dim(dist.X.test.train.U)[1]
  retval = apply(as.matrix(1:n.test),1,function(xx) {
             dist.max = sort(dist.X.test.train.L[xx,])[n.neigh]
             how.many.neigh=sum(dist.X.test.train.U[xx,]<=as.numeric(dist.max))
             return((how.many.neigh/n.neigh)*(n.L/n.U))
           })
  return(retval)
}

estimate_error_weights = function(pred.lab.test,pred.unlab.test,se=F)
{
  pred.lab.test[pred.lab.test<=0]     = 0
  pred.unlab.test[pred.unlab.test<=0] = 0

  output      = NULL
  output$mean = 1/2*mean(pred.lab.test^2)-mean(pred.unlab.test)
  if ( se==F ) return(output)

  # Standard Error Estimation
  n.L       = length(pred.lab.test)
  n.U       = length(pred.unlab.test)
  var.L     = var(pred.lab.test^2)
  var.U     = var(pred.unlab.test)
  output$se = sqrt(var.L/(4*n.L)+var.U/(n.U))
  return(output)
}


