
# Edit this to point to the R function file location
ABSOLUTE.PATH = "/Users/peterfreeman/PROJECTS/PHOTO_Z/CODE/"

source(paste(ABSOLUTE.PATH,"generic.R",sep=""))
source(paste(ABSOLUTE.PATH,"series.R",sep=""))
source(paste(ABSOLUTE.PATH,"kernn.R",sep=""))
source(paste(ABSOLUTE.PATH,"nn.R",sep=""))
source(paste(ABSOLUTE.PATH,"cde.R",sep=""))

# For creating figures akin to Figure 1 of Freeman, Izbicki, and Lee
plot_density = function(data.s,data.p,xlo,xhi,xlab,n=5000)
{
  density.spec=density(data.s, n = n)
  density.photo=density(data.p, n = n)
  maxY=max(c(density.spec$y,density.photo$y))
  plot(density.spec,col=2,lwd=3,xlim=c(xlo,xhi),ylim=c(0,maxY),xlab=xlab,main="",cex.main=1.7,cex.axis=1.5,cex.lab=1.8,cex=1.3)
  lines(density.photo,col=4,lwd=3,lty=2)
  legend("topright",c("Spec","Photo"),col=c(2,4),lwd=3,cex=1.5,bty="n",lty=1:2)
}

# Not explicitly necessary, but reduces the amount of repetition of code
# between entering spectroscopic and photometric data
input_data = function(fname,input_redshift=FALSE)
{
  redshift = NULL
  data     = read.table(fname,header=TRUE)
  cov      = cbind(
               data$u - data$g,
               data$g - data$r,
               data$r - data$i,
               data$i - data$z,
               data$z - data$y,
               data$r
             )
  if ( input_redshift == TRUE ) redshift = data$redshift
  return(list(cov=cov,redshift=redshift))
}

out = input_data("buzzard_spectroscopic.txt",input_redshift=TRUE)
cov.spec = out$cov
redshift = out$redshift
out = input_data("buzzard_photometric.txt",input_redshift=FALSE)
cov.photo = out$cov
rm(out)

#################
SHOW.PLOT = FALSE  # Set to true for the creation of some diagnostic plots.
#################

if ( SHOW.PLOT == TRUE ) { # The following plots could be created independently
  pdf("cde_figures.pdf")   # of the diagnostic plots...your call.
  plot_density(cov.spec$r,cov.photo$r,18,30,xlab="r magnitude")
  plot_density(cov.spec$u-cov.spec$g,cov.photo$u-cov.photo$g,-1,3.4,"u-g")
  plot_density(cov.spec$g-cov.spec$r,cov.photo$g-cov.photo$r,0,2.5,"g-r")
  plot_density(cov.spec$r-cov.spec$i,cov.photo$r-cov.photo$i,-0.2,1,"r-i")
  plot_density(cov.spec$i-cov.spec$z,cov.photo$i-cov.photo$z,-0.5,1,"i-z")
  plot_density(cov.spec$z-cov.spec$y,cov.photo$z-cov.photo$y,-0.5,1,"z-y")
}

############################################
out = compute_cde(cov.spec,cov.photo,redshift,n.L=15000,n.U=15000,min.U=1,
                  train.frac=7/15,val.frac=3/15,iseed=401,
                  COVARIATE.SHIFT=FALSE,SERIES=TRUE,KERNN=TRUE,NN=TRUE,
                  USE.Z.BOUNDS=TRUE,Z.MIN.SAMPLE=0,Z.MAX.SAMPLE=2,
                  NN.MIN=2,NN.MAX=30,NN.NUM=15,SHOW.PLOT=SHOW.PLOT,
              OUTPUT.DIR="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")
############################################

if ( SHOW.PLOT == TRUE ) dev.off()

