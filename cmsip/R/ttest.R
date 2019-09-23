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

write(sprintf('%d regions with depth >= %d\n', nrow(f), mindepth), stderr())

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
			write(rname, stderr())
			x=unname(unlist(tblA[rname, ]))
			y=unname(unlist(tblB[rname, ]))
			lfc=log2(
				(mean(x)+0.001) / (mean(y)+0.001)
				)
			if ((sd(x) + sd(y)) < 1e-10) {
				data.frame(
					symbol=rname
					, baseMean=mean(c(x, y))
					, lfc=lfc
					, statistic=mean(x) - mean(y)
					, pvalue=1e-300
					)
			} else {
				tmp=t.test(x, y, alternative='two.sided', var.equal=T)
				data.frame(
					symbol=rname
					, baseMean=mean(c(x, y))
					, lfc=lfc
					, statistic=tmp$statistic
					, pvalue=tmp$p.value
					)
			}
		}
		, mc.cores=numthreads
		)
	)
res$padj=p.adj(res$baseMean, res$pvalue)
res=na.omit(res)
res=res[order(res$padj),]
write.table(res, file=gzfile(outfile), quote=F, sep='\t', row.names=F)
