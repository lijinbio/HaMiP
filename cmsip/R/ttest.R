# vim: set noexpandtab tabstop=2:

f=read.table(infile, header=T, sep='\t', quote="", check.names=F, comment.char='#', stringsAsFactors=F)
rownames(f)=sprintf('%s:%d-%d', f$chrom, f$start, f$end)
f=f[, c(group1, group2), drop=F]
f=subset(f, rowSums(f)>=mindepth)

if(nrow(f)<1) {
	write(sprintf('Warning: no regions left under mindepth=%d\n', mindepth), stderr())
	q()
}

sf=read.table(sf_file, header=T, sep='\t', stringsAsFactors=F)[, c('sample_id', 'sizefactors'), drop=F]
rownames(sf) = sf[, 1]
sf[, 1]=NULL

f=as.data.frame(t(t(f) * sf[names(f), 1]))
f=split(f, rep(1:nsplit, ceiling(nrow(f)/nsplit))[1:nrow(f)])

res=do.call(rbind
	, parallel::mclapply(
		unname(f)
		, function (subf) {
			do.call(rbind
				, apply(
					subf
					, 1
					, function(rvalues) {
						x=rvalues[group1]
						y=rvalues[group2]
						lfc=log2(
							(mean(x)+0.001) / (mean(y)+0.001)
							)
						if ((sd(x) + sd(y)) < 1e-10) {
							data.frame(
								baseMean=mean(c(x, y))
								, lfc=lfc
								, statistic=mean(x) - mean(y)
								, pvalue=1e-300
								)
						} else {
							tmp=t.test(x, y, alternative='two.sided', var.equal=T)
							data.frame(
								baseMean=mean(c(x, y))
								, lfc=lfc
								, statistic=tmp$statistic
								, pvalue=tmp$p.value
								)
						}
					})
				)}
		, mc.cores=numthreads
		)
	)
res$padj=p.adj(res$baseMean, res$pvalue)
if (! keepNA) {
	res=na.omit(res)
}
res=res[order(res$padj),]
write.table(cbind(symbol=rownames(res), res), file=gzfile(outfile), quote=F, sep='\t', row.names=F)
