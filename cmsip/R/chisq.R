# vim: set noexpandtab tabstop=2:

f=read.table(infile, header=T, sep='\t', quote="", check.names=F, comment.char='#', stringsAsFactors=F)
symbls=sprintf('%s:%d-%d', f$chrom, f$start, f$end)
rownames(f)=symbls
f=f[, c(group1, group2), drop=F]
f=subset(f, rowSums(f)>=mindepth)
symbls=rownames(f)

if(nrow(f)<1) {
	write(sprintf('Warning: no regions left under mindepth=%d\n', mindepth), stderr())
	q()
}

sf=read.table(sf_file, header=T, sep='\t', stringsAsFactors=F)[, c('sample_id', 'sizefactors'), drop=F]
rownames(sf) = sf[, 1]
sf[, 1]=NULL

f=t(t(f) * sf[names(f), 1])
tblA=f[, group1, drop=F]
tblB=f[, group2, drop=F]

res=do.call(rbind
	, parallel::mclapply(
		symbls
		, function (rname) {
			x=unname(unlist(tblA[rname, ]))
			y=unname(unlist(tblB[rname, ]))
			lfc=log2(
				(mean(x)+0.001) / (mean(y)+0.001)
				)
			tmp=chisq.test(c(sum(x), sum(y))
				, p=c(length(x), length(y))
				, rescale.p=T)
			data.frame(
				symbol=rname
				, baseMean=mean(c(x, y))
				, lfc=lfc
				, statistic=tmp$statistic
				, pvalue=tmp$p.value
				)
		}
		, mc.cores=numthreads
		)
	)
res$padj=p.adj(res$baseMean, res$pvalue)
if (! keepNA) {
	res=na.omit(res)
}
res=res[order(res$padj),]
write.table(res, file=gzfile(outfile), quote=F, sep='\t', row.names=F)
