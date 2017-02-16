source("/Users/peterfreeman/PROJECTS/PHOTO_Z/CODE/generic.R")

combine_cde = function(filename.1,filename.2,
                       INPUT.DIR="./",OUTPUT.DIR=INPUT.DIR,
                       ALPHA.MIN=0,ALPHA.MAX=1,ALPHA.NUM=51)
{
  method = c(FALSE,FALSE,FALSE)
  load(paste(INPUT.DIR,filename.1,sep=""))
  if ( exists("series") ) {
    obj.1 = series; method[1] = TRUE
  } else if ( exists("kernn") ) {
    obj.1 = kernn; method[2] = TRUE
  } else if ( exists("nn") ) {
    obj.1 = nn; method[3] = TRUE
  }
  load(paste(INPUT.DIR,filename.2,sep=""))
  if ( exists("series") && method[1] == FALSE ) {
    obj.2 = series
  } else if ( exists("kernn") && method[2] == FALSE ) {
    obj.2 = kernn
  } else if ( exists("nn") && method[3] == FALSE ) {
    obj.2 = nn
  }

  alpha = seq(ALPHA.MIN,ALPHA.MAX,length.out=ALPHA.NUM)
  loss  = rep(NA,length(alpha))

  for ( ii in 1:length(alpha) ) {
    pred.val.U = alpha[ii]*obj.1$pred.val.U+(1-alpha[ii])*obj.2$pred.val.U
    pred.val.L = alpha[ii]*obj.1$pred.val.L+(1-alpha[ii])*obj.2$pred.val.L
    loss[ii]   = error_final_generic(pred.val.L,pred.val.U,
                                     z.val.L,weights.val.L,boot=F)$mean
  }

  best.alpha  = alpha[which.min(loss)]
  pred.test.L = best.alpha*obj.1$pred.test.L+(1-best.alpha)*obj.2$pred.test.L
  pred.test.U = best.alpha*obj.1$pred.test.U+(1-best.alpha)*obj.2$pred.test.U

  final.loss = error_final_generic(pred.test.L,pred.test.U,z.test.L,
                                   weights.test.L,boot=400)

  combined = list(pred.test.L=pred.test.L,pred.test.U=pred.test.U,
                  best.alpha=best.alpha,alpha=alpha,final.loss=final.loss)

  save(COVARIATE.SHIFT,combined,
       file=paste(OUTPUT.DIR,"cde_combined.Rdata",sep=""))
}

