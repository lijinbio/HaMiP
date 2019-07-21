# vim: set noexpandtab tabstop=2:

f=read.table(infile, header=T, sep='\t')
regions=paste(f$chrom, f$start, f$end, sep=':')
rownames(f)=regions

fg1=f[, group1]
fg2=f[, group2]

res=do.call(rbind
	, parallel::mclapply(
		regions
		, function(region) {
			x=unlist(fg1[region, ])
			y=unlist(fg2[region, ])
			lfc=log2(
				(mean(x)+0.001) / (mean(y)+0.001)
				)
			if ((sd(x) + sd(y)) == 0) {
				data.frame(
					lfc=lfc
					, statistic=mean(x) - mean(y)
					, pvalue=1e-300
					)
			} else {
				tmp=t.test(x, y, alternative='less', var.equal=T)
				data.frame(
					lfc=lfc
					, statistic=tmp$statistic
					, pvalue=tmp$p.value
					)
			}
		}
		, mc.cores=numthreads
		)
	)
res=cbind(chrom=f$chrom, start=f$start, end=f$end, fg1, fg2, res)
res$padj=p.adjust(res$pvalue, method='BH')
res=res[order(res$padj),]
write.table(res, file=outfile, quote=F, sep='\t', row.names=F)
