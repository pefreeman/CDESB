
coverage = function(prediction,redshift,plot=TRUE,
                    file.plot="coverage.pdf",output.dir="./")
{
  zGrid      = seq(0,1,length.out=dim(prediction)[2])
  zGrid.diff = as.matrix(dist(zGrid,method="manhattan"))
  prediction = t(apply(prediction,1,function(xx)xx/sum(xx))) # normalize
  prediction = t(apply(prediction,1,cumsum))
  m          = dim(prediction)[1]
  
  alpha = seq(0,1,length.out=101)

  alpha.hat = rep(0,101)
  
  ll = 1
  for ( ii in 1:m ) {
    if ( ii%%100 == 0 ) cat(ii," ")
    ll = ll+1
    diff = as.matrix(dist(prediction[ii,],method="manhattan")) 
    for ( jj in 2:100 ) {
      comp = diff>=alpha[jj]
      w = which(zGrid.diff==min(zGrid.diff[comp]),arr.ind=TRUE)
      if ( length(w) == 0 ) {
        alpha.hat[jj] = alpha.hat[jj] + 1
      } else {
        max.diff = 0
        max.lo = max.hi = -9
        for ( kk in 1:nrow(w) ) {
          if ( abs(prediction[ii,w[kk,1]] - 
                   prediction[ii,w[kk,2]]) > max.diff ) {
            max.diff = abs(prediction[ii,w[kk,1]] - prediction[ii,w[kk,2]])
            max.lo = min(w[kk,])
            max.hi = max(w[kk,])
          }
        }
        if ( zGrid[max.lo] <= redshift[ii] && zGrid[max.hi] >= redshift[ii] ) {
          alpha.hat[jj] = alpha.hat[jj] + 1
        }
      }
    }
  }
  alpha.hat = alpha.hat/m
  cat("\n")
  
  alpha.hat[101] = 1

  if ( plot == TRUE ) {
    pdf(paste(output.dir,file.plot,sep=""),height=6.5,width=6.5)
    plot(alpha,alpha.hat,typ="l",xlab="Expected Coverage",
         ylab="Observed Coverage",cex.axis=1.25,cex.lab=1.25,
         col=4,lwd=2,xlim=c(0,1),ylim=c(0,1))
    abline(a=0,b=1,lwd=2,col=2)
    dev.off()
  }
  return(list(alpha=alpha,alpha.hat=alpha.hat))
}

uniformity = function(prediction,redshift,test="chisq",plot=TRUE,
                      file.plot="uniformity.pdf",output.dir="./",
                      c=seq(0,1,length.out=101))
{
  prediction = t(apply(prediction,1,function(xx)xx/sum(xx))) # normalize
  prediction = t(apply(prediction,1,cumsum))
  m          = dim(prediction)[1]
  n          = dim(prediction)[2]
  num.bin    = length(c)-1
  c.hat      = rep(NA,num.bin)
  X          = rep(-9,m)
  for ( ii in 1:m ) {
    X[ii] = approx((0:(n-1))/(n-1),prediction[ii,],xout=redshift[ii])$y
  }
  if ( test == "chisq" || test == "qq" ) {
    X.bin = as.integer(floor(X*num.bin))+1
    w = which(X.bin>num.bin)
    X.bin[w] = num.bin
    X.bin.u = sort(unique(X.bin))
    obs = rep(0,num.bin)
    for ( ii in 1:num.bin ) {
      obs[ii] = length(which(X.bin==X.bin.u[ii]))
      if ( ii == 1 ) {
        c.hat[ii] = obs[ii]/m
      } else {
        c.hat[ii] = c.hat[ii-1]+obs[ii]/m
      }
    }
  }
  if ( plot == TRUE ) {
    pdf(paste(output.dir,file.plot,sep=""),height=6.5,width=6.5)
    if ( test == "chisq" || test == "qq" ) {
      xlab = "Expected Quantile"; ylab = "Observed Quantile"
      plot(c,c(0,c.hat),typ="l",xlab=xlab,ylab=ylab,cex.axis=1.25,cex.lab=1.25,
           col=4,lwd=2,xlim=c(0,1),ylim=c(0,1))
    } else {
      xlab = "Observed Percentile"; ylab = "Fraction of Observations"
      plot(ecdf(X),verticals=TRUE,col="blue",xlab=xlab,ylab=ylab,
           cex.lab=1.25,cex.axis=1.25,xlim=c(0,1),ylim=c(0,1),main="")
    }
    abline(a=0,b=1,lwd=2,col=2)
    dev.off()
  }
  chisq   = NA
  p.value = NA
  if ( test == "ks" ) {
    p.value = ks.test(X,"punif")$p.value
  } else if ( test == "ad" ) {
    p.value = ad.test(X)$p.value
  } else if ( test == "cvm" ) {
    p.value = cvm.test(X)$p.value
  } else if ( test == "chisq" || test == "qq" ) {
    chisq   = sum((obs - m/num.bin)^2/obs)
    p.value = 1-pchisq(chisq,num.bin)
  }
  return(list(p.value=p.value,chisq=chisq,c=c,c.hat=c.hat))
}
