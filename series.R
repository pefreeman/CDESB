source("/Users/peterfreeman/PROJECTS/PHOTO_Z/CODE/generic.R")

series = function(COVARIATE.SHIFT,z.train.L,z.val.L,z.test.L,
                  dist.train.L_train.L,dist.val.L_train.L,dist.val.U_train.L,
                  dist.test.L_train.L,dist.test.U_train.L,
                  weights.val.L,weights.test.L,
                  show.plot=TRUE,n.z.max=70,n.x.max=350,eps.min=0.1,
                  eps.max=100,eps.num=30,del.min=0,del.max=0.25,del.del=0.025,
                  z.min=0,z.max=1,z.num=201,boot.num=400)
{
  print("Series")

  eps.grid     = seq(log10(eps.min),log10(eps.max),length.out=eps.num)
  eps.grid     = exp(eps.grid)
  error        = rep(NA,length(eps.grid))

  for ( ii in 1:length(eps.grid) ) {
    cat(round(ii/length(eps.grid),2)," ")
    eps    = eps.grid[ii]
    obj.cden = cden_series(dist.train.L_train.L,z.train.L,
                           n.Z.max=n.z.max,n.X.max=n.x.max,
                           kernel.function=radial_kernel_distance,
                           extra.kernel=list("eps.val"=eps),
                           norm=NULL,system="Fourier")
    if ( COVARIATE.SHIFT == FALSE ) {
      obj.err = error_series(obj.cden,z.val.L,dist.val.L_train.L)
    } else {
      obj.err = error_series_cs(obj.cden,z.val.L,dist.val.L_train.L,
                                dist.val.U_train.L,weights.val.L)
    }
    error[ii] = obj.err$best.error
  }
  cat("\n")

  if ( show.plot == TRUE ) plot(eps.grid,error,pch=18)

  best.eps = eps.grid[which.min(error)]
  obj.cden = cden_series(dist.train.L_train.L,z.train.L,
                         n.Z.max=n.z.max,n.X.max=n.x.max,
                         kernel.function=radial_kernel_distance,
                         extra.kernel=list("eps.val"=best.eps),
                         norm=NULL,system="Fourier")
  if ( COVARIATE.SHIFT == FALSE ) {
    obj.err = error_series(obj.cden,z.val.L,dist.val.L_train.L)
  } else {
    obj.err = error_series_cs(obj.cden,z.val.L,dist.val.L_train.L,
                              dist.val.U_train.L,weights.val.L)
  }
  best.delta = choose_delta(obj.err,z.val.L,dist.val.L_train.L,
                            delta.grid=seq(del.min,del.max,del.del))

  pred.test.L = pred_density(
    obj.err,z.test.min=z.min,z.test.max=z.max,B=z.num,dist.test.L_train.L,
    prob.int=F,delta=best.delta$best.delta)

  pred.test.U = pred_density(
    obj.err,z.test.min=z.min,z.test.max=z.max,B=z.num,dist.test.U_train.L,
    prob.int=F,delta=best.delta$best.delta)

  pred.val.L  = pred_density(
    obj.err,z.test.min=z.min,z.test.max=z.max,B=z.num,dist.val.L_train.L,
    prob.int=F,delta=best.delta$best.delta)

  pred.val.U  = pred_density(
    obj.err,z.test.min=z.min,z.test.max=z.max,B=z.num,dist.val.U_train.L,
    prob.int=F,delta=best.delta$best.delta)

  final.loss  = error_final_generic(pred.test.L,pred.test.U,
                                    z.test.L,weights.test.L,
                                    boot=boot.num)

  return(list(pred.test.L=pred.test.L,pred.test.U=pred.test.U,
              pred.val.L=pred.val.L,pred.val.U=pred.val.U,
              final.loss=final.loss,object.error=obj.err,
              best.delta=best.delta$best.delta,best.eps=best.eps))
}

cden_series = function(distancesX,z,n.Z.max=20,n.X.max=NULL,
                       kernel.function=radialKernelDistance,
                       extra.kernel=list("eps.val"=1),
                       norm=NULL,system="Fourier")
{
  # estimates f(z|X)
  # z is assumed to be between 0 and 1
  # returns all coefficients up to n.Z.max, the maximum number of components
  # for Z, same for X
  kernel.matrix = kernel.function(distancesX,extra.kernel)
  if(any(is.na(kernel.matrix))) stop("Kernel with NA")
  n = length(distancesX[,1])

  norm.param                 = NULL
  norm.param$kernel.matrix.N = kernel.matrix
  kernel.matrix.N            = kernel.matrix
  if ( !is.null(norm) ) { # if data (X) needs to be normalized first
    if(norm=="symmetric") {
      norm.param      = symmetricNormalization(kernel.matrix)
      kernel.matrix.N = norm.param$kernel.matrix.N # normalized matrix
    } else {
      stop("Normalization rule not implemented!")
    }
  }

  #print("Computing Eigenvectors\n")

  if ( is.null(n.X.max) ) n.X.max = length(z)-10

  p            = 10
  Omega        = matrix(rnorm(n*(n.X.max+p),0,1),n,n.X.max+p)
  Z            = kernel.matrix.N%*%Omega
  Y            = kernel.matrix.N%*%Z
  Q            = qr(x=Y)
  Q            = qr.Q(Q)
  B            = t(Q)%*%Z%*%solve(t(Q)%*%Omega)
  eigenB       = eigen(B)
  lambda       = eigenB$values
  U            = Q%*%eigenB$vectors
  basis.X      = Re(sqrt(n)*U[,1:n.X.max])
  eigen.values = Re((lambda/n)[1:n.X.max])
  rm(Omega,Z,Y,Q,B,U,eigenB)
  gc()

  n.X = n.X.max

  #print("Done")

  # returns vector length(z) x n.Z with the basis for z.
  basis.Z = calculate_basis(z,n.Z.max,system)

  coeff = 1/n*t(t(basis.Z)%*%basis.X)

  obj                  = list()
  class(obj)           = "cDensity"
  obj$coeff            = matrix(coeff,n.X,n.Z.max)
  obj$system           = system
  obj$n.X              = n.X
  obj$n.Z              = n.Z.max
  obj$norm.param       = norm.param
  if ( is.null(norm) )   norm = 0
  obj$norm             = norm
  obj$eigen.X          = basis.X
  obj$eigen.values.X   = eigen.values
  obj$kernel.function  = kernel.function
  obj$extra.kernel     = extra.kernel
  return(obj)
}

error_series = function(obj,z.test,dist.test.train)
{
  # returns a n.X by n.Z matrix with the errors for the different possibilities
  # of number of components, plus the best combination
  if ( class(obj) != "cDensity" ) stop("Object should be of class cDensity")
 
  ker.new.old = obj$kernel.function(dist.test.train,obj$extra.kernel)
  if ( obj$norm != 0 ) {
    if ( obj$norm == "symmetric" ) {
      sqrt.col.means = obj$norm$sqrt.col.means
      sqrt.row.means = sqrt(rowMeans(ker.new.old))
      ker.new.old    = ker.new.old/(sqrt.row.means%*%t(sqrt.col.means))
    }
  }
  if ( any(is.na(ker.new.old)) ) stop("Kernel with NA")
 
  n.X = obj$n.X
  n.Z = obj$n.Z
  m   = dim(ker.new.old)[1] # New
  n   = dim(ker.new.old)[2] # Old
 
  # returns matrix length(z) x n.Z with the basis for z.
  basis.Z       = calculate_basis(z.test,n.Z,obj$system)
  eigen.vectors = obj$eigen.X
  eigen.values  = obj$eigen.values.X
  basis.X       = ker.new.old %*% eigen.vectors
  basis.X       = 1/n*basis.X*matrix(rep(1/eigen.values,m),m,n.X,byrow=T)
  rm(eigen.vectors)
 
  basis.psi.mean  = 1/m*t(t(basis.Z)%*%basis.X)
  W               = 1/m*t(basis.X)%*%basis.X

  prod.matrix     = lapply(as.matrix(1:n.Z),function(xx) {
                      aux.matrix=W*obj$coeff[,xx,drop=F]%*%t(obj$coeff[,xx,drop=F])
                      returnValue=diag(apply(t(apply(aux.matrix,1,cumsum)),2,cumsum))
                      return(returnValue[1:n.X])
                    })
  prod.matrix    = sapply(prod.matrix,function(xx)xx)
  D              = t(apply(prod.matrix,1,cumsum))
  rm(W)
  rm(prod.matrix)

  grid   = expand.grid(1:n.X,1:n.Z)
  errors = apply(grid,1,function(xx) {
             sBeta=1/2*D[xx[1],xx[2]]
             sLikeli=sum(obj$coeff[1:xx[1],1:xx[2]]*basis.psi.mean[1:xx[1],1:xx[2]])
             return(sBeta-sLikeli)}
           )
  errors = matrix(errors,n.X,n.Z)
  rm(D)
  rm(basis.psi.mean)

  point.min      = which(errors == min(errors,na.rm=T), arr.ind = TRUE)
  n.X.best       = (1:n.X)[point.min[1]]
  n.Z.best       = (1:n.Z)[point.min[2]]
  obj$n.X.best   = n.X.best
  obj$n.Z.best   = n.Z.best
  obj$errors     = errors
  obj$best.error = min(errors)
  return(obj)
}

error_series_cs = function(obj,z.val.L,dist.val.L_train.L,dist.val.U_train.L,
                           weight.z.val.L)
{
  # returns a n.X by n.Z matrix with the errors for the different possibilities
  # of number of components, plus the best combination
  if ( class(obj) != "cDensity" ) stop("Object should be of class cDensity")

  ker.test.L.train.L = obj$kernel.function(dist.val.L_train.L,obj$extra.kernel)
  ker.test.U.train.L = obj$kernel.function(dist.val.U_train.L,obj$extra.kernel)

  if ( any(is.na(ker.test.U.train.L)) ) stop("Kernel with NA")

  n.X       = obj$n.X
  n.Z       = obj$n.Z
  m.test.L  = dim(ker.test.L.train.L)[1] # New
  n.train.L = dim(ker.test.L.train.L)[2] # Old
  m.test.U  = dim(ker.test.U.train.L)[1] # New

  # returns matrix length(z)xn.Z with the basis for z.
  basis.Z   = calculate_basis(z.val.L,n.Z,obj$system)

  evecs = obj$eigen.X
  evals = obj$eigen.values.X

  basis.X.test.L = ker.test.L.train.L %*% evecs
  basis.X.test.L = 1/n.train.L*basis.X.test.L*
                     matrix(rep(1/evals,m.test.L),m.test.L,n.X,byrow=T)
  basis.X.test.U = ker.test.U.train.L %*% evecs
  basis.X.test.U = 1/n.train.L*basis.X.test.U*
                     matrix(rep(1/evals,m.test.U),m.test.U,n.X,byrow=T)
  rm(evecs)

  W = 1/m.test.U*t(basis.X.test.U)%*%basis.X.test.U

  matrix.beta    = matrix(weight.z.val.L,m.test.L,n.X)
  basis.psi.mean = 1/m.test.L*t(t(basis.Z)%*%(basis.X.test.L*(matrix.beta)))
  grid           = expand.grid(1:n.X,1:n.Z)

  prod.matrix = lapply(as.matrix(1:n.Z),function(xx) {
                  aux.matrix = W*obj$coeff[,xx,drop=F]%*%
                               t(obj$coeff[,xx,drop=F])
                  retval=diag(apply(t(apply(aux.matrix,1,cumsum)),2,cumsum))
                  return(retval[1:n.X])
                })
  prod.matrix = sapply(prod.matrix,function(xx)xx)
  D           = t(apply(prod.matrix,1,cumsum))
  rm(prod.matrix)

  errors = apply(grid,1,function(xx) {
             sBeta=1/2*D[xx[1],xx[2]]
             sLikeli=sum(obj$coeff[1:xx[1],1:xx[2]]*
                         basis.psi.mean[1:xx[1],1:xx[2]])
             return(sBeta-sLikeli)
           })
  errors = matrix(errors,n.X,n.Z)
  rm(D,W,basis.psi.mean)

  point.min      = which(errors == min(errors,na.rm=T), arr.ind = TRUE)
  obj$n.X.best   = (1:n.X)[point.min[1]]
  obj$n.Z.best   = (1:n.Z)[point.min[2]]
  obj$errors     = errors
  obj$best.error = min(errors)
  return(obj)
}

error_final = function(obj,z.test,dist.test.train,boot=F,delta=0)
{
  # estimates the loss of a given estimator, after fitting was done
  # boot==F if no error estimates, othewise boot is the number of
  #         bootstrap samples
  z.grid         = seq(0,1,0.005)
  pred           = pred_density(obj,z.test.min=0,z.test.max=1,B=length(z.grid),
                                dist.test.train,prob.int=F,delta=delta)
  col.means.comp = colMeans(pred^2)
  s.square       = mean(col.means.comp )
  n              = length(z.test)
  pred.obs       = apply(as.matrix(1:n),1,function(xx) {
                     index = which.min(abs(z.test[xx]-z.grid))
                     return(pred[xx,index])
                   })
  output         = NULL
  output$mean    = 1/2*s.square-mean(pred.obs)
  if ( boot==F ) return(output)

  # Bootstrap
  meanBoot = apply(as.matrix(1:boot),1,function(xx) {
               sample.boot    = sample(1:n,replace=T)
               pred.boot      = pred[sample.boot,]
               z.test.boot    = z.test[sample.boot]
               col.means.comp = colMeans(pred.boot^2)
               s.square       = mean(col.means.comp )
               pred.obs       = apply(as.matrix(1:n),1,function(xx) {
                                  index=which.min(abs(z.test.boot[xx]-z.grid))
                                  return(pred.boot[xx,index])
                                })
               return(1/2*s.square-mean(pred.obs))
             })
  output$seBoot = sqrt(var(meanBoot))
  return(output)
}

calculate_basis = function(z,n.Z,system)
{
  if ( system == "cosine" ) {
    basis.Z = apply(as.matrix(1:(n.Z-1)),1,function(xx)sqrt(2)*cos(xx*pi*z))
    basis.Z = cbind(rep(1,length(z)),basis.Z)
    return(basis.Z)
  } else if ( system == "Fourier" ) {
    sin.basis.Z = apply(as.matrix(1:round((n.Z)/2)),1,
                        function(xx) sqrt(2)*sin(2*xx*pi*z)
                       )
    cos.basis.Z = apply(as.matrix(1:round((n.Z)/2)),1,
                        function(xx) sqrt(2)*cos(2*xx*pi*z)
                       )
    basis.Z = matrix(NA,length(z),2*round((n.Z)/2))
    basis.Z[,seq(from=1,length.out=dim(sin.basis.Z)[2],by=2)] = sin.basis.Z
    basis.Z[,seq(from=2,length.out=dim(cos.basis.Z)[2],by=2)] = cos.basis.Z
    basis.Z = cbind(rep(1,length(z)),basis.Z)
    basis.Z = basis.Z[,1:n.Z]
    return(basis.Z)
  } else {
    stop("Error: basis unsupported.")
  }
}

choose_delta = function(obj,z.val,dist.val.train,delta.grid=seq(0,0.3,0.04)) {
  error = errorSE = rep(NA,length(delta.grid))
  for ( ii in 1:length(delta.grid) ) {
    cat(round(ii/length(delta.grid),2)," ")
    estimate.errors = error_final(obj,z.val,
                      dist.val.train,boot=100,delta=delta.grid[ii])
    error[ii]   = estimate.errors$mean
    errorSE[ii] = estimate.errors$seBoot
  }
  cat("\n")
  plot(delta.grid,error+errorSE,type="l",lwd=2,
       ylim=c(min(error),max(error+errorSE)))
  lines(delta.grid,error,lwd=3)
  whichMin              = (1:length(error))[error==min(error)]
  whichMin              = max(whichMin)
  best.delta            = delta.grid[whichMin]
  best.delta.1SD        = max(delta.grid[error<=(error+errorSE)[whichMin]])
  retval                = NULL
  retval$best.delta     = best.delta
  retval$best.delta.1SD = best.delta.1SD
  return(retval)
}

pred_density = function(obj,z.test.min=0,z.test.max=1,B=1000,
                        dist.test.train,prob.int=F,confidence=0.95,delta=0)
{
  # predict density at points z.test=seq(z.test.min,z.test.max,length.out=B)
  # and xTest; delta is the treshhold to remove bumps
  #print("Densities normalized to integrate 1 in the range of z given!")

  if ( class(obj) != "cDensity" ) stop("Object should be of class cDensity")
  z.grid = seq(from=z.test.min,to=z.test.max,length.out=B)

  ker.new.old = obj$kernel.function(dist.test.train,obj$extra.kernel)
  if ( obj$norm != 0 ) {
    if ( obj$norm == "symmetric" ) {
      sqrt.col.means = obj$norm.param$sqrt.col.means
      sqrt.row.means = sqrt(rowMeans(ker.new.old))
      ker.new.old    = ker.new.old/(sqrt.row.means%*%t(sqrt.col.means))
    }
  }
  n.X.best = obj$n.X
  n.Z.best = obj$n.Z
  if ( !is.null(obj$n.X.best) ) { # if CV was done
    n.X.best = obj$n.X.best
    n.Z.best = obj$n.Z.best
  }
  m = dim(ker.new.old)[1] # New
  n = dim(ker.new.old)[2] # Old

  # returns matrix length(z) x n.Z with the basis for z.
  basis.Z = calculate_basis(z.grid,n.Z.best,obj$system)

  eigen.vectors = obj$eigen.X[,1:n.X.best,drop=F]
  eigen.values  = obj$eigen.values.X[1:n.X.best]
  basis.X       = ker.new.old %*% eigen.vectors
  basis.X       = 1/n*basis.X*matrix(rep(1/eigen.values,m),m,n.X.best,byrow=T)
  n.test.Z      = length(z.grid) # how many test points for z
  n.test.X      = dim(ker.new.old)[1] # how many test points for x
  grid          = expand.grid(1:n.test.X,1:n.test.Z)
  coeff         = obj$coeff[1:n.X.best,1:n.Z.best,drop=F]
  estimates     = t(apply(basis.X,1,function(yy) {
                      apply(basis.Z,1,function(xx)
                        sum(outer(yy,xx)*coeff)
                      )
                    }))
  bin.size      = (z.test.max-z.test.min)/(B+1)
  estimates     = t(apply(estimates,1,function(xx)
                      normalize_density(bin.size,xx,delta)
                    ))
  if ( !prob.int ) return(estimates)

  # Gives threshold on density corresponding to probability interval
  thresh         = apply(estimates,1,function(xx)
                     findThreshold(bin.size,xx,confidence)
                   )
  obj            = list()
  obj$estimates  = estimates
  obj$thresh.int = thresh
  return(obj)
}

radial_kernel = function(distances,extra.kernel=list("eps.val"=1))
{
  # Given the distances and the bandwidth eps.val,
  # computes the matrix of radial kernel
  return(exp(-distances^2/(4*extra.kernel$eps.val)))
}

radial_kernel_distance=function(distances,extra.kernel=list("eps.val"=1))
{
  # Given the distances and the bandwidth eps.val,
  # computes the matrix of radial kernel
  return(exp(-distances^2/(4*extra.kernel$eps.val)))
}

