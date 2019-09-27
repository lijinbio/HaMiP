# From DESeq2::pvalueAdjustment()
p.adj=function(
	stats
	, pvalue
	, alpha=0.1
	, method='BH'
	) {
	theta=seq(0, 0.99, length=50)
	filtPadj=genefilter::filtered_p(
		filter=stats
		, test=pvalue
		, theta=theta
		, method=method
		) 
	numRej=colSums(filtPadj<alpha, na.rm=T)
	lo.fit=lowess(numRej ~ theta, f=1/5)
	if (max(numRej) <= 10) {
		j=1
	} else { 
		residual=if (all(numRej==0)) { 0 } else { numRej[numRej > 0] - lo.fit$y[numRej > 0] }
		thresh=max(lo.fit$y) - sqrt(mean(residual^2))
		j=if (any(numRej > thresh)) { which(numRej > thresh)[1] } else { 1 }
	}
	filtPadj[, j]
}
