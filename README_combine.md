
R functions for combining conditional density estimates based on 
"A unified framework for constructing, tuning and assessing photometric
redshift density estimates in a selection bias setting" by
Freeman, Izbicki and Lee (2017).

---

To combine two CDEs, edit and use "fe_combine.R" (meaning "front end to the 
combine conditional density estimations code"). (Note that the code is currently
limited to the case where two CDEs are combined.) The example that is provided 
assumes that we are combining the outputs from using the Series and kerNN 
estimators. 

The arguments to combine_cde are:

filename.1	some pair combination of cde_series.Rdata, cde_kernn.Rdata,
filename.2	and cde_nn.Rdata
INPUT.DIR	where the cde_*.Rdata files exist
OUTPUT.DIR	where to place cde_combined.Rdata
ALPHA.MIN	vector of linear combination coefficients to try when
ALPHA.MAX	determining the optimal combination
ALPHA.NUM

The output file cde_combined.Rdata contains the following:

COVARIATE.SHIFT	a record of whether this was TRUE or FALSE
combined	a list similar to the series, kernn, and nn lists
		documented in README_cde.md...in particular, it contains
		the matrix pred.test.L, each row of which is p(z) for a
		labeled test set galaxy


