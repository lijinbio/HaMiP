#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.1.2"

import os
import sys
import argparse
import yaml
import subprocess

def runcmd(cmd, log=subprocess.PIPE, echo=False):
	if echo:
		print('Running: ' + cmd)
	try:
		cp=subprocess.run('bash -c "%s"' % cmd, universal_newlines=True, shell=True, stdout=log, stderr=subprocess.STDOUT)
		if cp.returncode != 0:
			print('Error: %s failed.' % cmd, vars(cp), sep='\n', file=sys.stderr)
			sys.exit(-1)
	except OSError as e:
		print("Execution failed: ", e, file=sys.stderr)
		sys.exit(-1)
	return cp

def runcmdsh(cmd, log=subprocess.PIPE, echo=False):
	if echo:
		print('Running: ' + cmd)
	try:
		cp=subprocess.run(cmd, universal_newlines=True, shell=True, stdout=log, stderr=subprocess.STDOUT)
		if cp.returncode != 0:
			print('Error: %s failed.' % cmd, vars(cp), sep='\n', file=sys.stderr)
			sys.exit(-1)
	except OSError as e:
		print("Execution failed: ", e, file=sys.stderr)
		sys.exit(-1)
	return cp

def bsmap_runcmd(fname, refenece, numthread, outfile, verbose=False):
	if os.path.exists(outfile):
		return
	runcmd('mkdir -p ' + os.path.dirname(outfile), echo=verbose)
	cmd = 'bsmap' + \
		' -a ' + fname + \
		' -d ' + refenece + \
		' -R -n 1 -r 0 ' + \
		' -p ' + str(numthread) + \
		' -o ' + outfile
	runcmd(cmd, log=open(outfile+".stdout", 'w+'), echo=verbose)

def bsmap_ref(config, reference):
	outbasedir=os.path.join(config['resultdir'], 'bsmap', reference)
	for sampleinfo in config['sampleinfo']:
		outfile=os.path.join(outbasedir, sampleinfo['sampleid'] + '.bam')
		if os.path.exists(outfile):
			continue
		if len(sampleinfo['filenames']) > 1:
			files = ''
			for f in sampleinfo['filenames']:
				bname=os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0];
				fname=os.path.join(config['aligninfo']['fastqdir'], f)
				singlefile=os.path.join(outbasedir, 'single', bname + '.bam')
				bsmap_runcmd(fname, config['aligninfo'][reference], config['aligninfo']['numthreads'], singlefile)
				files += ' ' + singlefile
			runcmd('samtools merge ' + outfile + ' ' + files, echo=config['aligninfo']['verbose'])
		else:
			fname=os.path.join(config['aligninfo']['fastqdir'], sampleinfo['filenames'][0])
			bsmap_runcmd(fname, config['aligninfo'][reference], config['aligninfo']['numthreads'], outfile, config['aligninfo']['verbose'])

import re
def bsmap_stat_parse(infile):
	with open(infile) as f:
		dstr=f.read()
	totalreads=int(re.search('total reads: (\d+)', dstr).groups()[0])
	alignedreads=int(re.search('aligned reads: (\d+)', dstr).groups()[0])
	uniquereads=int(re.search('unique reads: (\d+)', dstr).groups()[0])
	return (totalreads, alignedreads, uniquereads)

def bsmap_stat(config, reference):
	basedir=os.path.join(config['resultdir'], 'bsmap')
	stats = {}
	for sampleinfo in config['sampleinfo']:
		if len(sampleinfo['filenames']) > 1:
			totalr=0
			alignedr=0
			uniquer=0
			for fname in sampleinfo['filenames']:
				bname=os.path.splitext(os.path.splitext(os.path.basename(fname))[0])[0];
				f=os.path.join(basedir, reference, 'single', bname + '.bam.stdout')
				ftotalr, falignedr, funiquer = bsmap_stat_parse(f)
				totalr += ftotalr
				alignedr += falignedr
				uniquer += funiquer
			stats[sampleinfo['sampleid']] = (totalr, alignedr, uniquer)
		else:
			f=os.path.join(basedir, reference, sampleinfo['sampleid'] + '.bam.stdout')
			stats[sampleinfo['sampleid']] = bsmap_stat_parse(f)
	return stats

def bsmap(config):
	if config['aligninfo']['verbose']:
		print('==>bsmap<==')
	bsmap_ref(config, 'reference')
	bsmap_ref(config, 'spikein')
	mpstat = {}
	mpstat['reference'] = bsmap_stat(config, 'reference')
	mpstat['spikein'] = bsmap_stat(config, 'spikein')
	if config['aligninfo']['verbose']:
		print(mpstat)
	return mpstat

def removeCommonReads_runcmd(infile1, infile2, outfile1, outfile2, verbose=False):
	bin=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'perl', 'removeCommonRead.pl')
	runcmd('mkdir -p ' + os.path.dirname(outfile1) + ' ' + os.path.dirname(outfile2), echo=verbose)
	cmd = bin + ' ' + infile1 + ' ' + infile2 + \
		' >(' + 'samtools view -bS - ' + ' -o ' + outfile1 + ' 2>/dev/null)' \
		' >(' + 'samtools view -bS - ' + ' -o ' + outfile2 + ' 2>/dev/null)'
	cp=runcmd(cmd, echo=verbose)
	return int(cp.stdout.strip())

def removeCommonReads(config):
	if config['aligninfo']['verbose']:
		print('==>removeCommonReads<==')
	inbasedir=os.path.join(config['resultdir'], 'bsmap')
	outbasedir=os.path.join(config['resultdir'], 'removeCommonReads')
	comm={}
	for sampleinfo in config['sampleinfo']:
		refinfile=os.path.join(inbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkinfile=os.path.join(inbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		refoutfile=os.path.join(outbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkoutfile=os.path.join(outbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		comm[sampleinfo['sampleid']] = removeCommonReads_runcmd(refinfile, spkinfile, refoutfile, spkoutfile, config['aligninfo']['verbose'])
	if config['aligninfo']['verbose']:
		print(comm)
	return comm

def totalwigsums_n(f):
	cmd = "bedtools genomecov -ibam %s -bg | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { sum=0 } { sum += $4*($3-$2) } END { print sum }'" % f
	cp=runcmdsh(cmd)
	return int(cp.stdout.strip())

def totalwigsums(config):
	if config['aligninfo']['verbose']:
		print('==>totalwigsums<==')
	inbasedir=os.path.join(config['resultdir'], 'removeCommonReads')
	twss = {}
	for reference in ['reference', 'spikein']:
		tws = {}
		for sampleinfo in config['sampleinfo']:
			f=os.path.join(inbasedir, reference, sampleinfo['sampleid'] + '.bam')
			tws[sampleinfo['sampleid']] = totalwigsums_n(f)
		twss[reference] = tws
	if config['aligninfo']['verbose']:
		print(twss)
	return twss

import statistics
def estimateSizeFactors(tws, verbose=False):
	if verbose:
		print('==>estimateSizeFactors<==')
	mws = statistics.median(tws.values())
	sizefactors = {id : mws / ws for id, ws in tws.items()}
	if verbose:
		print(sizefactors)
	return sizefactors

def normalizetwsref(tws, sizefactors, verbose=False):
	if verbose:
		print('==>normalizetwsref<==')
	twsn = {id: ws * sizefactors[id] for id, ws in tws.items()}
	if verbose:
		print(twsn)
	return twsn

def saveQCstats(config, statfile, qcstats):
	with open(statfile, 'w+') as f:
		print('\t'.join((
			'sample_id'
			, 'total'
			, 'unique_ref'
			, 'ref/total'
			, 'unique_spk'
			, 'spk/total'
			, 'comm'
			, 'comm/total'
			, 'comm/unique_ref'
			, 'twss_spk'
			, 'sizefactors'
			, 'twss_ref'
			, 'twss_ref_norm'
			)), file=f)
		for sampleinfo in config['sampleinfo']:
			sampleid = sampleinfo['sampleid']
			total = qcstats['mpstat']['reference'][sampleid][0]
			unique_ref = qcstats['mpstat']['reference'][sampleid][2]
			unique_spk = qcstats['mpstat']['spikein'][sampleid][2]
			comm = qcstats['comm'][sampleid]
			twss_spk = qcstats['twss']['spikein'][sampleid]
			sizefactors = qcstats['sizefactors'][sampleid]
			twss_ref = qcstats['twss']['reference'][sampleid]
			twss_ref_norm = qcstats['twsrefnorm'][sampleid]
			print('\t'.join(map(str
				, (sampleid
					, total
					, unique_ref
					, '{:.2%}'.format(unique_ref/total)
					, unique_spk
					, '{:.2%}'.format(unique_spk/total)
					, comm
					, '{:.2%}'.format(comm/total)
					, '{:.2%}'.format(comm/unique_ref)
					, twss_spk
					, '{:.2f}'.format(sizefactors)
					, twss_ref
					, '{:.0f}'.format(twss_ref_norm)
					)
				)), file=f)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
def barplot(config, tws):
	outfile=os.path.join(config['resultdir'], config['aligninfo']['barplotinfo']['outfile'])
	plt.figure(figsize=(config['aligninfo']['barplotinfo']['width'], config['aligninfo']['barplotinfo']['height']))
	ax=plt.axes()
	ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: '{:,}'.format(int(x))))
	plt.bar(*zip(*tws.items()), width=0.5, color='blue')
	plt.title('Normalized total wigsum')
	runcmd('mkdir -p ' + os.path.dirname(outfile))
	plt.savefig(outfile, bbox_inches='tight')

def removedupref(config):
	if config['aligninfo']['verbose']:
		print('==>removedupref<==')
	indir=os.path.join(config['resultdir'], 'removeCommonReads', 'reference')
	outdir=os.path.join(config['resultdir'], 'removedupref')
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(indir, sampleinfo['sampleid'] + '.bam')
			outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bam')
			runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['aligninfo']['verbose'])
			runcmd('samtools rmdup %s %s' % (infile, outfile), echo=config['aligninfo']['verbose'])

def bamtobed(config):
	if config['genomescaninfo']['verbose']:
		print('==>bamtobed<==')
	indir=os.path.join(config['resultdir'], 'removedupref')
	outdir=os.path.join(config['resultdir'], 'bamtobed')
	for sampleinfo in config['sampleinfo']:
		infile=os.path.join(indir, sampleinfo['sampleid'] + '.bam')
		outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bed')
		if os.path.exists(outfile):
			continue
		runcmd('mkdir -p ' + os.path.dirname(outfile))
		runcmd('bedtools bamtobed -i %s > %s' % (infile, outfile))
	return outdir

def readextension(config):
	if config['genomescaninfo']['verbose']:
		print('==>readextension<==')
	indir=os.path.join(config['resultdir'], 'bamtobed')
	outdir=os.path.join(config['resultdir'], 'readextension_frag'+str(config['genomescaninfo']['fragsize']))
	for sampleinfo in config['sampleinfo']:
		infile=os.path.join(indir, sampleinfo['sampleid'] + '.bed')
		outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bed')
		if os.path.exists(outfile):
			continue
		runcmd('mkdir -p ' + os.path.dirname(outfile))
		cmd="awk -v FS='\\t' -v OFS='\\t' -v fragsize=%s -e '{ if ($6==\"+\") { $3=$2+fragsize } else if ($6==\"-\") { $2=$3-fragsize; if($2<0) { $2=0 } } print }' < %s > %s" % (
				config['genomescaninfo']['fragsize'], infile, outfile)
		runcmdsh(cmd, config['genomescaninfo']['verbose'])
	return outdir

def fetchChromSizes(config):
	outfile=os.path.join(config['resultdir'], config['genomescaninfo']['referencename'] + '_genomefile.txt')
	if os.path.exists(outfile):
		return outfile
	runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['genomescaninfo']['verbose'])
	runcmd('fetchChromSizes %s > %s' % (config['genomescaninfo']['referencename'], outfile), echo=config['genomescaninfo']['verbose'])
	return outfile

def makewindows(genomefile, windowsize, windowfile):
	runcmd('mkdir -p ' + os.path.dirname(windowfile))
	runcmd('bedtools makewindows -g %s -w %s | sort -k 1,1 -k 2,2n -k 3,3n > %s' % (genomefile, windowsize, windowfile))

import tempfile
def tabulatereadcounts(config, windowfile, beddir, counttablefile):
	if config['genomescaninfo']['verbose']:
		print('==>tabulatereadcounts<==')
	cntdir=tempfile.TemporaryDirectory(dir=config['resultdir']).name
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(beddir, sampleinfo['sampleid'] + '.bed')
			outfile=os.path.join(cntdir, sampleinfo['sampleid'] + '.bedgraph')
			runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['genomescaninfo']['verbose'])
			cmd = "bedtools coverage -a %s -b %s -counts | awk -v FS='\\t' -v OFS='\\t' -e '$4>0' > %s" % (windowfile, infile, outfile)
			runcmdsh(cmd, config['genomescaninfo']['verbose'])
	sampleids=[sampleinfo['sampleid'] for sampleinfo in config['sampleinfo']]
	fs=[os.path.join(cntdir, id+'.bedgraph') for id in sampleids]
	runcmd('mkdir -p ' + os.path.dirname(counttablefile), echo=config['genomescaninfo']['verbose'])
	runcmd("bedtools unionbedg -i %s -header -names %s | gzip -n > %s" % (' '.join(fs), ' '.join(sampleids), counttablefile), echo=config['genomescaninfo']['verbose'])
	for sampleinfo in config['sampleinfo']:
			runcmd('rm -f ' + os.path.join(cntdir, sampleinfo['sampleid'] + '.bedgraph'), echo=config['genomescaninfo']['verbose'])
	runcmd('rmdir --ignore-fail-on-non-empty ' + cntdir, echo=config['genomescaninfo']['verbose'])

def tabulatemeanwig(config, windowfile, genomefile, beddir, counttablefile):
	if config['genomescaninfo']['verbose']:
		print('==>tabulatemeanwig<==')
	cntdir=tempfile.TemporaryDirectory(dir=config['resultdir']).name
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(beddir, sampleinfo['sampleid'] + '.bed')
			covfile=os.path.join(cntdir, sampleinfo['sampleid'] + '.genomecov.bedgraph')
			runcmd('mkdir -p ' + os.path.dirname(covfile), echo=config['genomescaninfo']['verbose'])
			runcmd("bedtools genomecov -i %s -g %s -bg > %s" % (infile, genomefile, covfile), echo=config['genomescaninfo']['verbose'])
			outfile=os.path.join(cntdir, sampleinfo['sampleid'] + '.bedgraph')
			runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['genomescaninfo']['verbose'])
			cmd = "bedtools intersect -a %s -b %s -wo | awk -v FS='\\t' -v OFS='\\t' -e '{ print $1, $2, $3, $7*$8/%d }' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools groupby -g 1,2,3 -c 4 -o sum > %s" % (windowfile, covfile, config['genomescaninfo']['windowsize'], outfile)
			runcmdsh(cmd, config['genomescaninfo']['verbose'])
	sampleids=[sampleinfo['sampleid'] for sampleinfo in config['sampleinfo']]
	fs=[os.path.join(cntdir, id+'.bedgraph') for id in sampleids]
	runcmd('mkdir -p ' + os.path.dirname(counttablefile), echo=config['genomescaninfo']['verbose'])
	runcmd("bedtools unionbedg -i %s -header -names %s | gzip -n > %s" % (' '.join(fs), ' '.join(sampleids), counttablefile), echo=config['genomescaninfo']['verbose'])
	for sampleinfo in config['sampleinfo']:
			runcmd('rm -f ' + os.path.join(cntdir, sampleinfo['sampleid'] + '.genomecov.bedgraph'), echo=config['genomescaninfo']['verbose'])
			runcmd('rm -f ' + os.path.join(cntdir, sampleinfo['sampleid'] + '.bedgraph'), echo=config['genomescaninfo']['verbose'])
	runcmd('rmdir --ignore-fail-on-non-empty ' + cntdir, echo=config['genomescaninfo']['verbose'])

def swapdict(d):
	nd = {}
	for k, v in d.items():
		if v in nd:
			nd[v].append(k)
		else:
			nd[v]=[k]
	return nd

def ttest(config, statfile, counttablefile, testfile):
	if config['dhmrinfo']['verbose']:
		print('==>t.test<==')
	sampleid2group={sampleinfo['sampleid']:sampleinfo['group'] for sampleinfo in config['sampleinfo']}
	group2sampleid=swapdict(sampleid2group)
	group1=group2sampleid[config['groupinfo']['group1']]
	group2=group2sampleid[config['groupinfo']['group2']]
	g1str = 'c(' + ', '.join(["'" + name + "'" for name in group1]) + ')'
	g2str = 'c(' + ', '.join(["'" + name + "'" for name in group2]) + ')'
	adjscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'p.adj.R')
	rscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'ttest.R')
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"numthreads=%s\" -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\" -e \"source('%s')\"" % (
			config['dhmrinfo']['numthreads'], counttablefile, statfile, config['dhmrinfo']['mindepth'], g1str, g2str, testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', adjscript, rscript
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

def chisq(config, statfile, counttablefile, testfile):
	if config['dhmrinfo']['verbose']:
		print('==>chisq.test<==')
	sampleid2group={sampleinfo['sampleid']:sampleinfo['group'] for sampleinfo in config['sampleinfo']}
	group2sampleid=swapdict(sampleid2group)
	group1=group2sampleid[config['groupinfo']['group1']]
	group2=group2sampleid[config['groupinfo']['group2']]
	g1str = 'c(' + ', '.join(["'" + name + "'" for name in group1]) + ')'
	g2str = 'c(' + ', '.join(["'" + name + "'" for name in group2]) + ')'
	adjscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'p.adj.R')
	rscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'chisq.R')
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"numthreads=%s\" -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\" -e \"source('%s')\"" % (
			config['dhmrinfo']['numthreads'], counttablefile, statfile, config['dhmrinfo']['mindepth'], g1str, g2str, testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', adjscript, rscript
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

def gtest(config, statfile, counttablefile, testfile):
	if config['dhmrinfo']['verbose']:
		print('==>gtest<==')
	sampleid2group={sampleinfo['sampleid']:sampleinfo['group'] for sampleinfo in config['sampleinfo']}
	group2sampleid=swapdict(sampleid2group)
	group1=group2sampleid[config['groupinfo']['group1']]
	group2=group2sampleid[config['groupinfo']['group2']]
	g1str = 'c(' + ', '.join(["'" + name + "'" for name in group1]) + ')'
	g2str = 'c(' + ', '.join(["'" + name + "'" for name in group2]) + ')'
	adjscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'p.adj.R')
	rscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'gtest.R')
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"numthreads=%s\" -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\" -e \"source('%s')\"" % (
			config['dhmrinfo']['numthreads'], counttablefile, statfile, config['dhmrinfo']['mindepth'], g1str, g2str, testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', adjscript, rscript
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

def nbtest(config, statfile, counttablefile, testfile):
	if config['dhmrinfo']['verbose']:
		print('==>nbtest<==')
	sampleid2group={sampleinfo['sampleid']:sampleinfo['group'] for sampleinfo in config['sampleinfo']}
	group2sampleid=swapdict(sampleid2group)
	group1=group2sampleid[config['groupinfo']['group1']]
	group2=group2sampleid[config['groupinfo']['group2']]
	g1str = 'c(' + ', '.join(["'" + name + "'" for name in group1]) + ')'
	g2str = 'c(' + ', '.join(["'" + name + "'" for name in group2]) + ')'
	rscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'nbtest.R')
	windowsize=0
	if not config['genomescaninfo']['readscount']:
		windowsize=config['genomescaninfo']['windowsize']
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"windowsize=%s\" -e \"condA='%s'\" -e \"condB='%s'\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\"" % (
			counttablefile, statfile, config['dhmrinfo']['mindepth'], g1str, g2str, windowsize, config['groupinfo']['group1'], config['groupinfo']['group2'], testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', rscript
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

def align_run(config):
	statfile = os.path.join(config['resultdir'], 'qcstats.txt')
	if 'statfile' in config['aligninfo']:
		statfile = config['aligninfo']['statfile']
	if not os.path.exists(statfile):
		qcstats = {}
		qcstats['mpstat'] = bsmap(config)
		qcstats['comm'] = removeCommonReads(config)
		qcstats['twss'] = totalwigsums(config)
		qcstats['sizefactors'] = estimateSizeFactors(qcstats['twss']['spikein'], config['aligninfo']['verbose'])
		qcstats['twsrefnorm'] = normalizetwsref(qcstats['twss']['reference'], qcstats['sizefactors'], config['aligninfo']['verbose'])
		saveQCstats(config, statfile, qcstats)
		barplot(config, qcstats['twsrefnorm'])
		removedupref(config)
	return statfile

def genomescan_run(config):
	counttablefile=os.path.join(config['resultdir'], 'counttable.txt.gz')
	if 'counttablefile' in config['genomescaninfo']:
		counttablefile = config['genomescaninfo']['counttablefile']
	if not os.path.exists(counttablefile):
		beddir=bamtobed(config)
		if config['genomescaninfo']['readextension']:
			beddir=readextension(config)
		genomefile=fetchChromSizes(config)
		windowfile=config['genomescaninfo']['windowfile']
		if not os.path.exists(windowfile):
			makewindows(genomefile, config['genomescaninfo']['windowsize'], windowfile)
		if config['genomescaninfo']['readscount']:
			tabulatereadcounts(config, windowfile, beddir, counttablefile)
		else:
			tabulatemeanwig(config, windowfile, genomefile, beddir, counttablefile)
	return counttablefile

def mergedhmr(config, testfile, outfile):
	runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['dhmrinfo']['verbose'])
	hyperfile=outfile+'.hyper.bed'
	hypofile=outfile+'.hypo.bed'

	qcol=7
	if config['dhmrinfo']['method'] in ['ttest', 'chisq', 'gtest']:
		qcol=6

	cmd = "(printf '%%s\\n' chrom start end $(zcat %s | head -n 1 | cut -f 2-) | paste -s -d $'\\t'; zcat %s | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { getline } $%s<%s && $3>0 { n=split($1, a, \"[:-]\"); for (i=1; i<=n; i++) { printf a[i] OFS } for (j=2; j<=NF; j++) { printf $j ((j<NF)?OFS:ORS) } }' | sort -k 1,1 -k 2,2n -k 3,3n) > %s" % (
				testfile, testfile, qcol, config['dhmrinfo']['qthr'], hyperfile
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])
	cmd = "(printf '%%s\\n' chrom start end $(zcat %s | head -n 1 | cut -f 2-) | paste -s -d $'\\t'; zcat %s | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { getline } $%s<%s && $3<0 { n=split($1, a, \"[:-]\"); for (i=1; i<=n; i++) { printf a[i] OFS } for (j=2; j<=NF; j++) { printf $j ((j<NF)?OFS:ORS) } }' | sort -k 1,1 -k 2,2n -k 3,3n) > %s" % (
				testfile, testfile, qcol, config['dhmrinfo']['qthr'], hypofile
			)

	if config['dhmrinfo']['method'] in ['ttest', 'chisq', 'gtest']:
		runcmd("(head -n 1 %s; (bedtools merge -i %s -d %s -c 4,5,6,7,8 -o max,max,max,min,min -header | tail -n +2; bedtools merge -i %s -d %s -c 4,5,6,7,8 -o max,min,max,min,min -header | tail -n +2) | sort -k 8,8g ) | gzip -n > %s" % (
			hyperfile, hyperfile, config['dhmrinfo']['maxdistance'], hypofile, config['dhmrinfo']['maxdistance'], outfile
			), echo=config['dhmrinfo']['verbose'])
	else:
		runcmd("bedtools merge -i %s -d %s -c 4,5,6,7,8,9 -o max,absmax,absmax,absmax,min,min -header | gzip -n > %s" % (tmpfile, config['dhmrinfo']['maxdistance'], outfile), echo=config['dhmrinfo']['verbose'])
		runcmd("(head -n 1 %s; (bedtools merge -i %s -d %s -c 4,5,6,7,8,9 -o max,max,max,max,min,min -header | tail -n +2; bedtools merge -i %s -d %s -c 4,5,6,7,8,9 -o max,min,max,min,min,min -header | tail -n +2) | sort -k 9,9g ) | gzip -n > %s" % (
			hyperfile, hyperfile, config['dhmrinfo']['maxdistance'], hypofile, config['dhmrinfo']['maxdistance'], outfile
			), echo=config['dhmrinfo']['verbose'])
	runcmd('rm -f %s %s' % (hyperfile, hypofile), echo=config['dhmrinfo']['verbose'])

def dhmr_run(config, statfile, counttablefile):
	testfile=os.path.join(config['resultdir'], 'testfile.txt.gz')
	if 'testfile' in config['dhmrinfo']:
		testfile = config['dhmrinfo']['testfile']
	if not os.path.exists(testfile):
		runcmd('mkdir -p ' + os.path.dirname(testfile), echo=config['dhmrinfo']['verbose'])
		if config['dhmrinfo']['method'] == 'ttest':
			ttest(config, statfile, counttablefile, testfile)
		elif config['dhmrinfo']['method'] == 'chisq':
			chisq(config, statfile, counttablefile, testfile)
		elif config['dhmrinfo']['method'] == 'gtest':
			gtest(config, statfile, counttablefile, testfile)
		else:
			nbtest(config, statfile, counttablefile, testfile)
	dhmrfile=os.path.join(config['resultdir'], 'dhmr.txt.gz')
	if 'dhmrfile' in config['dhmrinfo']:
		dhmrfile = config['dhmrinfo']['dhmrfile']
	if not os.path.exists(dhmrfile):
		mergedhmr(config, testfile, dhmrfile)

def run(config):
	statfile=align_run(config)
	counttablefile=genomescan_run(config)
	dhmr_run(config, statfile, counttablefile)

def updatedefs(data, pars):
	obj=data
	path=pars[0].split('.')
	for k in path[:-1]:
		if k not in obj:
			obj[k] = {}
		obj=obj[k]
	v=pars[1]
	if v.replace('.' , '', 1).isdigit():
		if '.' in v:
			v=float(v)
		else:
			v=int(v)
	elif v in ['True', 'T', 'true', 't']:
		v=True
	elif v in ['False', 'F', 'false', 'f']:
		v=False
	obj[path[-1]]=v

import re
import pprint
def main():
	parser = argparse.ArgumentParser(
		description='CMS-IP sequencing analysis'
		)
	parser.add_argument('-c', '--config'
		, type=argparse.FileType('r')
		, required=True
		, help='Configuration file in YAML format.'
		)
	parser.add_argument('-D'
		, action='append'
		, nargs='+'
		, type=lambda kv: re.split('=', kv)
		, help='Define variable=value to suppress configuration file. e.g.\n"-D dhmrinfo.verbose=False"'
		)
	parser.add_argument('-v', '--version'
		, action='version'
		, version='%(prog)s ' + __version__
		)
	args = parser.parse_args()
	config=yaml.load(args.config, Loader=yaml.FullLoader)
	if args.D is not None:
		for vars in args.D:
			for vs in vars:
				updatedefs(config, vs)
	pprint.pprint(config)
	run(config)

if __name__ == "__main__":
	main()
