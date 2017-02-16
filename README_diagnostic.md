
R functions for performing diagnostic tests of conditional density estimates,
based on "A unified framework for constructing, tuning and assessing 
photometric redshift density estimates in a selection bias setting" by
Freeman, Izbicki and Lee (2017).

---

To perform diagnostic tests, edit and use "fe_diagnostic.R" (meaning "front
end to the diagnostic tests code"). The example that is provided assumes
that we wish to perform diagnostic tests of the combined Series-kerNN 
estimator, dubbed "combined" in the input Rdata file.

There are two diagnostic functions:

1) uniformity

Tests for uniformity of values of hat(F)(z_true|x) for each galaxy.

The arguments are:

prediction	matrix of p(z) estimates, rows=galaxies and cols=redshifts
redshift	true redshift corresponding to each row of prediction
test		one of the following: "ks", "ad", "cvm", "chisq"
		"ks"	tests whether percentiles hat(F)(z_true|x) for all 
			galaxies are sampled from a uniform distribution, 
			via the KS test
			returns: p.value
		"ad"	same as "ks," but applies Anderson-Darling test
		"cvm"	same as "ks," but applies Cramer-von Mises test
		"chisq"	bins the observed percentiles hat(F)(z_true|x), then 
		or "qq"	tests whether the number of observed percentiles in 
			each bin is consistent with uniformity via the 
			chi-square statistic
			returns: chisq, p.value, vector of expected 
				 quantiles c, vector of observed quantiles 
				 c.hat
		******	=> NOTE: THIS OPTION PRODUCES A QQ PLOT if plot==TRUE
plot		if TRUE, generates a pdf plot with name file.plot in output.dir
file.plot
output.dir
c		sequence of quantiles for binning hat(F)(z|x); only used
		for test=="chisq"

The returned values are listed above, under "test."

2) coverage

Tests for coverage: how often do intervals containing a fraction alpha of the
integrated hat(f)(z|x) contain the true redshift? For instance, do intervals
constructed assuming alpha = 0.8 contain the true redshift 80% of the time?

The arguments are:

prediction      matrix of p(z) estimates, rows=galaxies and cols=redshifts
redshift        true redshift corresponding to each row of prediction
plot            if TRUE, generates a pdf plot with name file.plot in output.dir
file.plot
output.dir

The returned values:

alpha		sequence of alpha values
alpha.hat	sequence of fractions of redshifts lying within the
		constructed intervals with integrated probability alpha

