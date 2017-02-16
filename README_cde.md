
R functions for creating conditional density estimates based on 
"A unified framework for constructing, tuning and assessing photometric
redshift density estimates in a selection bias setting" by
Freeman, Izbicki and Lee (2017).

---

To construct CDEs, edit and use "fe_cde.R" (meaning "front end to the main
conditional density estimation code"). The example that is provided is
based on having two ASCII files, one with labeled data (i.e., covariate
measurements plus redshifts) and one with unlabeled data (i.e., covariate
measurements only). One would have to replace this code entirely given other
data formats (e.g., FITS or HDF5). The necessary variables to set values for
when inputting data are 

cov.spec 	matrix of labeled data covariates 
		(rows--each galaxy, columns--each covariate measurement)
cov.photo	same as cov.spec, but for unlabeled data
redshift	vector of redshifts (the labels) for labeled data

One then calls compute_cde:

compute_cde = function(cov.spec,cov.photo,redshift,n.L=15000,n.U=15000,min.U=1,
                       train.frac=7/15,val.frac=3/15,iseed=401,
                       COVARIATE.SHIFT=FALSE,SERIES=FALSE,KERNN=FALSE,NN=FALSE,
                       USE.Z.BOUNDS=TRUE,Z.MIN.SAMPLE=0,Z.MAX.SAMPLE=2,
                       NN.MIN=2,NN.MAX=30,NN.NUM=15,SHOW.PLOT=SHOW.PLOT,
               OUTPUT.DIR="/Users/peterfreeman/PROJECTS/PHOTO_Z/DC1/BUZZARD/")

The optional arguments are as follows:

n.L		the number of galaxies to extract from cov.spec (or cov.photo)
n.U		for use in analyses (because there are matrix computations 
		involved, one cannot generally use all the data in cov.spec and
		cov.photo...for typical desktop computers, the upper limit 
 		is ~ 10,000)
min.U           when resampling, this is the minimum number of unlabeled data
                that must be among the k nearest neighbors of a considered
                labeled point
train.frac	the fraction of the n.L (and n.U) galaxies that are to be used
val.frac 	for training and validation. The remaining galaxies are used
		as the test set
iseed		random-number generator seed

COVARIATE.SHIFT	TRUE if importance weights are to be computed and risk functions
		based on covariate shift are to be evaluated 
SERIES		TRUE if CDEs based on the Series method are to be computed
KERNN		TRUE if CDEs based on the kerNN method are to be computed
NN		TRUE if CDEs based on the NN method are to be computed
USE.Z.BOUNDS	TRUE if you wish to set the range of redshifts via 
Z.MIN.SAMPLE	Z.MIN.SAMPLE and Z.MAX.SAMPLE; if FALSE, redshifts are mapped
Z.MAX.SAMPLE	to the interval [0,1] via
		z.map = (z-min(z))/(max(z)-min(z))
NN.MIN		when determining importance weights, evaluate the risk function
NN.MAX		for each number of nearest neighbors in the sequence
NN.NUM		seq(NN.MIN,NN.MAX,length.out=NN.NUM)
SHOW.PLOT	TRUE if you wish to output diagnostic plots; if you wish the 
		plots to go to a specific plotting device (e.g., a PDF file), 
                that device must be specified within fe_cde.R *before* the
		call to compute_cde
OUTPUT.DIR	where to place the output files:
		if SERIES == TRUE : cde_series.Rdata
		if KERNN  == TRUE : cde_kernn.Rdata
		if NN     == TRUE : cde_nn.Rdata

Variables that are recorded in all Rdata output files are:

COVARIATE.SHIFT TRUE is covariate shift mitigation is applied
z.train.L	redshifts for the labeled data used for training, validation,
z.val.L		and testing...these lie in the interval [0,1]
z.test.L
z.train.L.orig	same as above, but the original, untransformed values
z.val.L.orig
z.test.L.orig
weights.train.L	importance weight estimates for the labeled training, validation,
weights.val.L	and test data, as well as the unlabeled test data
weights.test.L
weights.test.U
cov.train.L	the covariates for the labeled training set
sample.L	the rows of cov.spec and cov.photo used as the labeled
sample.U	and unlabeled datasets
ran.perm.L	the random ordering of the integers 1 -> sample.L (or sample.U)
ran.perm.U	used when assigning the labeled data to training, validation,
		and test sets...to determine which rows of the original input
		data correspond to these sets, use, e.g,	
			train.L.indices = sample.L[ran.perm.L[1:n.train]]
			val.L.indices   = sample.L[ran.perm.L[(n.train+1):
				  			      (n.train+n.val]]
			test.L.indices  = sample.L[ran.perm.L[(n.train+n.val+1):
							       n.L]]
n.train.L	the numbers of galaxies assigned to the training, validation, and
n.val.L		test sets for the labeled and unlabeled data
n.test.L
n.train.U
n.val.U
n.test.U
file		output file name

Variables that differ from Rdata file to Rdata file:

series		a list containing the following elements:
 		pred.val.L	matrices containing p(z) on each row
 		pred.val.U	(rows--galaxies,columns--by default z=[0,1]
		pred.test.L	in (Delta)z = 0.005 increments)
		pred.test.U
		final.loss	the evaluated loss function; compare this
				number to those from kernn, nn, etc.
		other outputs	kept to generate predictions for the
				full cov.photo file via predict.R

kernn		a list containing the following elements:
                same as for series, except note that the value of best.nn should
		be compared against the default upper limit of 30; if 
		kernn$best.nn == 30, edit cde.R to increase the
		upper limit and rerun the analysis

nn		a list containing the following elements:
		same as for series, except note that the value for best.nn should
		be compared against the default upper limit of 90; if
		nn$best.nn == 90, edit cde.R to increase the
                upper limit and rerun the analysis

