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

if (windowsize>0) {
	f=round(f*windowsize)
}

sample_table=data.frame(
	sample_id=c(group1, group2)
	, condition=c(
		rep(condA, length(group1))
		, rep(condB, length(group2))
		)
	)
design=as.formula(sprintf('~%s', names(sample_table)[[2]]))
suppressPackageStartupMessages(library(DESeq2))
dds=DESeqDataSetFromMatrix(
	countData=f
	, colData=sample_table
	, design=design
	)

sf=read.table(sf_file, header=T, sep='\t', stringsAsFactors=F)[, c('sample_id', 'sizefactors'), drop=F]
rownames(sf) = sf[, 1]
sf[, 1]=NULL

sizeFactors(dds) <- 1.0 / sf[names(f), 1]
if (windowsize>0) {
	sizeFactors(dds) <- 1.0 * windowsize / sf[names(f), 1]
}
res=data.frame(
	results(
		DESeq(dds)
		, contrast=c(names(sample_table)[[2]], condA, condB)
		)
	)
if (! keepNA) {
	res=na.omit(res)
}
res=res[order(res$padj),]
write.table(cbind(symbol=rownames(res), res), file=gzfile(outfile), quote = F, sep = '\t', row.names = F)
