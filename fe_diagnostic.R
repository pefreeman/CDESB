library(goftest)  # install this package if not installed already...

# Edit the following to point to your installation of diagnostic.R
source("/Users/peterfreeman/PROJECTS/PHOTO_Z/CODE/diagnostic.R")

# Here I assume that the files I need are in the current working directory.
load("cde_kernn.Rdata")		# gives the z.test.L
load("cde_combined.Rdata")	# gives combined$pred.test.L

predict  = combined$pred.test.L
redshift = z.test.L

# Coverage Plot -- as currently coded, slow
out.cov = coverage(predict,redshift,
          output.dir="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")

# Tests of Uniformity -- as currently coded, faster
out.unif = uniformity(predict,redshift,test="ks",file.plot="ks.pdf",
           output.dir="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")
cat("The p value for the KS test is: ",out.unif$p.value,"\n")

out.unif = uniformity(predict,redshift,test="ad",file.plot="ad.pdf",
           output.dir="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")
cat("The p value for the AD test is: ",out.unif$p.value,"\n")

out.unif = uniformity(predict,redshift,test="cvm",file.plot="cvm.pdf",
           output.dir="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")
cat("The p value for the CvM test is: ",out.unif$p.value,"\n")

out.unif = uniformity(predict,redshift,test="qq",file.plot="qq.pdf",
           output.dir="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")
cat("The chi-square value is:         ",out.unif$chisq,"\n")
cat("The p value for the GoF test is: ",out.unif$p.value,"\n")

