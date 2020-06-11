#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.2.8"

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

def makedirectory(path, echo=False):
	runcmd('mkdir -p %s' % (os.path.dirname(path) or '.'), echo=echo)

def lnfile(infile, outfile, verbose=False):
	if os.path.exists(outfile):
		return
	makedirectory(outfile, verbose)
	runcmd('ln -sf %s %s' % (infile, outfile), echo=verbose)

def bsmap_runcmd(fname, reference, numthread, outfile, verbose=False):
	if os.path.exists(outfile):
		return
	makedirectory(outfile, verbose)
	cmd = 'bsmap -a %s -d %s -n 1 -r 0 -S 1234 -p %d -o %s' % (fname, reference, numthread, outfile)
	runcmd(cmd, log=open(outfile+".stdout", 'w+'), echo=verbose)

def bsmap_runcmd_pe(fname1, fname2, reference, numthread, outfile, verbose=False):
	if os.path.exists(outfile):
		return
	makedirectory(outfile, verbose)
	cmd = 'bsmap -a %s -b %s -d %s -n 1 -r 0 -S 1234 -p %d -o %s' % (fname1, fname2, reference, numthread, outfile)
	runcmd(cmd, log=open(outfile+".stdout", 'w+'), echo=verbose)

def bsmap_ref(config, reference):
	outbasedir=os.path.join(config['resultdir'], 'bsmap', reference)
	for sampleinfo in config['sampleinfo']:
		outfile=os.path.join(outbasedir, sampleinfo['sampleid'] + '.bam')
		if os.path.exists(outfile):
			continue
		if 'inputbam' in config['aligninfo'] and config['aligninfo']: # input BAM files
			if len(sampleinfo[reference]) > 1:
				runcmd('samtools merge ' + outfile + ' ' + ' '.join(sampleinfo[reference]), echo=config['aligninfo']['verbose'])
			else:
				fname=sampleinfo[reference][0]
				lnfile(fname, outfile, config['aligninfo']['verbose'])
		else:
			if len(sampleinfo['filenames']) > 1:
				files = ''
				for f in sampleinfo['filenames']:
					bname=os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0];
					fname=f
					singlefile=os.path.join(outbasedir, 'single', bname + '.bam')
					bsmap_runcmd(fname, config['aligninfo'][reference], config['aligninfo']['numthreads'], singlefile)
					files += ' ' + singlefile
				runcmd('samtools merge ' + outfile + ' ' + files, echo=config['aligninfo']['verbose'])
			else:
				fname=sampleinfo['filenames'][0]
				bsmap_runcmd(fname, config['aligninfo'][reference], config['aligninfo']['numthreads'], outfile, config['aligninfo']['verbose'])

import re
def bsmap_stat_parse(infile):
	if not os.path.exists(infile):
		return (0, 0, 0)
	with open(infile) as f:
		dstr=f.read()
	totalreads=int(re.search('total reads: (\d+)', dstr).groups()[0])
	alignedreads=int(re.search('aligned reads: (\d+)', dstr).groups()[0])
	uniquereads=int(re.search('unique reads: (\d+)', dstr).groups()[0])
	return (totalreads, alignedreads, uniquereads)

def bam_numreads(infile):
	if not os.path.exists(infile):
		return 0
	cmd='samtools view -c %s' % infile
	cp=runcmdsh(cmd)
	return int(cp.stdout.strip())

def bsmap_stat(config, reference):
	basedir=os.path.join(config['resultdir'], 'bsmap')
	stats = {}
	if 'inputbam' in config['aligninfo'] and config['aligninfo']:
		for sampleinfo in config['sampleinfo']:
			if 'filenames' in sampleinfo and len(sampleinfo['filenames']) > 1:
				totalr=0
				for fname in sampleinfo['filenames']:
					bname=os.path.splitext(os.path.splitext(os.path.basename(fname))[0])[0];
					f=os.path.join(basedir, reference, 'single', bname + '.bam')
					totalr += bam_numreads(f)
				stats[sampleinfo['sampleid']] = (totalr, totalr, totalr)
			else:
				f=os.path.join(basedir, reference, sampleinfo['sampleid'] + '.bam')
				totalr=bam_numreads(f)
				stats[sampleinfo['sampleid']] = (totalr, totalr, totalr)
	else:
		for sampleinfo in config['sampleinfo']:
			if 'filenames' in sampleinfo and len(sampleinfo['filenames']) > 1:
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
	mpstat = {}
	bsmap_ref(config, 'reference')
	mpstat['reference'] = bsmap_stat(config, 'reference')
	if config['aligninfo']['usespikein']:
		bsmap_ref(config, 'spikein')
		mpstat['spikein'] = bsmap_stat(config, 'spikein')
	if config['aligninfo']['verbose']:
		print(mpstat)
	return mpstat

def mcall_stat_parse(infile):
	if not os.path.exists(infile):
		return 0
	with open(infile) as f:
		dstr=f.read()
	return float(re.search('bisulfite conversion ratio = ([\d.]+)', dstr).groups()[0])

def mcall_runcmd(infile, outdir, sampleid, reference, numthread, verbose=False):
	linkfile=os.path.join(outdir, sampleid + '.bam')
	statfile=linkfile+'_stat.txt'
	if os.path.exists(statfile):
		return mcall_stat_parse(statfile)
	lnfile(os.path.abspath(infile), linkfile)
	cmd = 'cd %s && mcall -m %s -r %s --sampleName %s -p %s' % (os.path.dirname(linkfile), os.path.basename(linkfile), reference, sampleid, numthread)
	runcmd(cmd, log=open(linkfile+".stdout", 'w+'), echo=verbose)
	return mcall_stat_parse(statfile)

def mcall_ref(config, reference):
	inbasedir=os.path.join(config['resultdir'], 'bsmap', reference)
	outbasedir=os.path.join(config['resultdir'], 'mcall', reference)
	stats = {}
	for sampleinfo in config['sampleinfo']:
		infile=os.path.join(inbasedir, sampleinfo['sampleid'] + '.bam')
		outdir=os.path.join(outbasedir, sampleinfo['sampleid'])
		stats[sampleinfo['sampleid']] = mcall_runcmd(infile
				, outdir
				, sampleinfo['sampleid']
				, config['aligninfo'][reference]
				, config['aligninfo']['numthreads']
				, config['aligninfo']['verbose']
				)
	return stats

def mcall(config):
	if config['aligninfo']['verbose']:
		print('==>mcall<==')
	mcstat = {}
	mcstat['reference'] = mcall_ref(config, 'reference')
	if config['aligninfo']['usespikein']:
		mcstat['spikein'] = mcall_ref(config, 'spikein')
	if config['aligninfo']['verbose']:
		print(mcstat)
	return mcstat

def removeCommonReads_runcmd(infile1, infile2, outfile1, outfile2, verbose=False):
	bin=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'perl', 'removeCommonRead.pl')
	makedirectory(outfile1)
	makedirectory(outfile2)
	outsam1=outfile1+".sam"
	outsam2=outfile2+".sam"
	cmd = bin + ' ' + infile1 + ' ' + infile2 + ' ' + outsam1 + ' ' + outsam2
	cp=runcmd(cmd, echo=verbose)
	runcmd('samtools view -bS -o %s %s' % (outfile1, outsam1))
	runcmd('samtools view -bS -o %s %s' % (outfile2, outsam2))
	runcmd('rm -f %s %s' % (outsam1, outsam2))
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
	inbasedir=os.path.join(config['resultdir'], 'removeCommonReads' if config['aligninfo']['usespikein'] else 'bsmap')
	twss = {}
	for reference in ['reference', 'spikein'] if config['aligninfo']['usespikein'] else ['reference']:
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
		if config['aligninfo']['usespikein']:
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
				, 'bcr_ref'
				, 'bcr_spk'
				)), file=f)
			for sampleinfo in config['sampleinfo']:
				sampleid = sampleinfo['sampleid']
				total = qcstats['mpstat']['reference'][sampleid][0] if sampleid in qcstats['mpstat']['reference'] else 0
				unique_ref = qcstats['mpstat']['reference'][sampleid][2] if sampleid in qcstats['mpstat']['reference'] else 0
				unique_spk = qcstats['mpstat']['spikein'][sampleid][2] if sampleid in qcstats['mpstat']['reference'] else 0
				comm = qcstats['comm'][sampleid]
				twss_spk = qcstats['twss']['spikein'][sampleid]
				sizefactors = qcstats['sizefactors'][sampleid]
				twss_ref = qcstats['twss']['reference'][sampleid]
				twss_ref_norm = qcstats['twsrefnorm'][sampleid]
				bcr_ref = qcstats['mcstat']['reference'][sampleid]
				bcr_spk = qcstats['mcstat']['spikein'][sampleid]
				print('\t'.join(map(str
					, (sampleid
						, total
						, unique_ref
						, '{:.2%}'.format(unique_ref/total) if total>0 else 'NA'
						, unique_spk
						, '{:.2%}'.format(unique_spk/total) if total>0 else 'NA'
						, comm
						, '{:.2%}'.format(comm/total) if total>0 else 'NA'
						, '{:.2%}'.format(comm/unique_ref) if unique_ref>0 else 'NA'
						, twss_spk
						, '{:.2f}'.format(sizefactors)
						, twss_ref
						, '{:.0f}'.format(twss_ref_norm)
						, '{:.6f}'.format(bcr_ref)
						, '{:.6f}'.format(bcr_spk)
						)
					)), file=f)
		else:
			print('\t'.join((
				'sample_id'
				, 'total'
				, 'unique_ref'
				, 'ref/total'
				, 'twss_ref'
				, 'sizefactors'
				, 'twss_ref_norm'
				, 'bcr_ref'
				)), file=f)
			for sampleinfo in config['sampleinfo']:
				sampleid = sampleinfo['sampleid']
				total = qcstats['mpstat']['reference'][sampleid][0] if sampleid in qcstats['mpstat']['reference'] else 0
				unique_ref = qcstats['mpstat']['reference'][sampleid][2] if sampleid in qcstats['mpstat']['reference'] else 0
				twss_ref = qcstats['twss']['reference'][sampleid]
				sizefactors = qcstats['sizefactors'][sampleid]
				twss_ref_norm = qcstats['twsrefnorm'][sampleid]
				bcr_ref = qcstats['mcstat']['reference'][sampleid]
				print('\t'.join(map(str
					, (sampleid
						, total
						, unique_ref
						, '{:.2%}'.format(unique_ref/total) if total>0 else 'NA'
						, twss_ref
						, '{:.2f}'.format(sizefactors)
						, '{:.0f}'.format(twss_ref_norm)
						, '{:.6f}'.format(bcr_ref)
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
	makedirectory(outfile)
	plt.savefig(outfile, bbox_inches='tight')

def removedupref(config):
	if config['aligninfo']['verbose']:
		print('==>removedupref<==')
	indir=os.path.join(config['resultdir'], 'removeCommonReads' if config['aligninfo']['usespikein'] else 'bsmap', 'reference')
	outdir=os.path.join(config['resultdir'], 'removedupref')
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(indir, sampleinfo['sampleid'] + '.bam')
			outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bam')
			makedirectory(outfile, config['aligninfo']['verbose'])
			runcmd('samtools rmdup -s %s %s' % (infile, outfile), echo=config['aligninfo']['verbose'])

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
		makedirectory(outfile)
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
		makedirectory(outfile)
		cmd="awk -v FS='\\t' -v OFS='\\t' -v fragsize=%s -e '{ if ($6==\"+\") { $3=$2+fragsize } else if ($6==\"-\") { $2=$3-fragsize; if($2<0) { $2=0 } } print }' < %s > %s" % (
				config['genomescaninfo']['fragsize'], infile, outfile)
		runcmdsh(cmd, config['genomescaninfo']['verbose'])
	return outdir

def fetchChromSizes(config):
	outfile=os.path.join(config['resultdir'], config['genomescaninfo']['referencename'] + '_genomefile.txt')
	if os.path.exists(outfile):
		return outfile
	makedirectory(outfile, config['genomescaninfo']['verbose'])
	runcmd('fetchChromSizes %s > %s' % (config['genomescaninfo']['referencename'], outfile), echo=config['genomescaninfo']['verbose'])
	return outfile

def makewindows(genomefile, windowsize, windowfile):
	makedirectory(windowfile)
	runcmd('bedtools makewindows -g %s -w %s | sort -k 1,1 -k 2,2n -k 3,3n > %s' % (genomefile, windowsize, windowfile))

import tempfile
def tabulatereadcounts(config, windowfile, beddir, counttablefile):
	if config['genomescaninfo']['verbose']:
		print('==>tabulatereadcounts<==')
	cntdir=tempfile.TemporaryDirectory(dir=config['resultdir']).name
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(beddir, sampleinfo['sampleid'] + '.bed')
			outfile=os.path.join(cntdir, sampleinfo['sampleid'] + '.bedgraph')
			makedirectory(outfile, config['genomescaninfo']['verbose'])
			cmd = "bedtools coverage -a %s -b %s -counts | awk -v FS='\\t' -v OFS='\\t' -e '$4>0' > %s" % (windowfile, infile, outfile)
			runcmdsh(cmd, config['genomescaninfo']['verbose'])
	sampleids=[sampleinfo['sampleid'] for sampleinfo in config['sampleinfo']]
	fs=[os.path.join(cntdir, id+'.bedgraph') for id in sampleids]
	makedirectory(counttablefile, config['genomescaninfo']['verbose'])
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
			makedirectory(covfile, config['genomescaninfo']['verbose'])
			runcmd("bedtools genomecov -i %s -g %s -bg > %s" % (infile, genomefile, covfile), echo=config['genomescaninfo']['verbose'])
			outfile=os.path.join(cntdir, sampleinfo['sampleid'] + '.bedgraph')
			makedirectory(outfile, config['genomescaninfo']['verbose'])
			cmd = "bedtools intersect -a %s -b %s -wo | awk -v FS='\\t' -v OFS='\\t' -e '{ print $1, $2, $3, $7*$8/%d }' | sort -k 1,1 -k 2,2n -k 3,3n | bedtools groupby -g 1,2,3 -c 4 -o sum > %s" % (windowfile, covfile, config['genomescaninfo']['windowsize'], outfile)
			runcmdsh(cmd, config['genomescaninfo']['verbose'])
	sampleids=[sampleinfo['sampleid'] for sampleinfo in config['sampleinfo']]
	fs=[os.path.join(cntdir, id+'.bedgraph') for id in sampleids]
	makedirectory(counttablefile, config['genomescaninfo']['verbose'])
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
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"numthreads=%s\" -e \"nsplit=%s\" -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\" -e \"source('%s')\"" % (
			config['dhmrinfo']['numthreads'], config['dhmrinfo']['nsplit'], counttablefile, statfile, config['dhmrinfo']['meandepth']*(len(group1)+len(group2)), g1str, g2str, testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', adjscript, rscript
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
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"numthreads=%s\" -e \"nsplit=%s\" -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\" -e \"source('%s')\"" % (
			config['dhmrinfo']['numthreads'], config['dhmrinfo']['nsplit'], counttablefile, statfile, config['dhmrinfo']['meandepth']*(len(group1)+len(group2)), g1str, g2str, testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', adjscript, rscript
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
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"numthreads=%s\" -e \"nsplit=%s\" -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\" -e \"source('%s')\"" % (
			config['dhmrinfo']['numthreads'], config['dhmrinfo']['nsplit'], counttablefile, statfile, config['dhmrinfo']['meandepth']*(len(group1)+len(group2)), g1str, g2str, testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', adjscript, rscript
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
	if 'readscount' in config['genomescaninfo'] and not config['genomescaninfo']['readscount']:
		windowsize=config['genomescaninfo']['windowsize']
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"infile='%s'\" -e \"sf_file='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"windowsize=%s\" -e \"condA='%s'\" -e \"condB='%s'\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\"" % (
			counttablefile, statfile, config['dhmrinfo']['meandepth']*(len(group1)+len(group2)), g1str, g2str, windowsize, config['groupinfo']['group1'], config['groupinfo']['group2'], testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', rscript
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

def nbtest_sf(config, counttablefile, testfile):
	if config['dhmrinfo']['verbose']:
		print('==>nbtest_sf<==')
	sampleid2group={sampleinfo['sampleid']:sampleinfo['group'] for sampleinfo in config['sampleinfo']}
	group2sampleid=swapdict(sampleid2group)
	group1=group2sampleid[config['groupinfo']['group1']]
	group2=group2sampleid[config['groupinfo']['group2']]
	g1str = 'c(' + ', '.join(["'" + name + "'" for name in group1]) + ')'
	g2str = 'c(' + ', '.join(["'" + name + "'" for name in group2]) + ')'
	rscript=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'R', 'nbtest_sf.R')
	windowsize=0
	if 'readscount' in config['genomescaninfo'] and not config['genomescaninfo']['readscount']:
		windowsize=config['genomescaninfo']['windowsize']
	cmd = "R --slave --no-save --no-restore --no-init-file -e \"infile='%s'\" -e \"mindepth=%d\" -e \"group1=%s\" -e \"group2=%s\" -e \"windowsize=%s\" -e \"condA='%s'\" -e \"condB='%s'\" -e \"outfile='%s'\" -e \"keepNA=%s\" -e \"source('%s')\"" % (
			counttablefile, config['dhmrinfo']['meandepth']*(len(group1)+len(group2)), g1str, g2str, windowsize, config['groupinfo']['group1'], config['groupinfo']['group2'], testfile, 'T' if config['dhmrinfo']['keepNA'] else 'F', rscript
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

def align_run(config):
	statfile = config['aligninfo']['statfile'] if 'statfile' in config['aligninfo'] else os.path.join(config['resultdir'], 'qcstats.txt')
	if not os.path.exists(statfile):
		qcstats = {}
		qcstats['mpstat'] = bsmap(config)
		qcstats['mcstat'] = mcall(config)
		if config['aligninfo']['usespikein']:
			qcstats['comm'] = removeCommonReads(config)
		qcstats['twss'] = totalwigsums(config)
		qcstats['sizefactors'] = estimateSizeFactors(
				qcstats['twss']['spikein'] if config['aligninfo']['usespikein'] else qcstats['twss']['reference']
				, config['aligninfo']['verbose'])
		qcstats['twsrefnorm'] = normalizetwsref(qcstats['twss']['reference'], qcstats['sizefactors'], config['aligninfo']['verbose'])
		saveQCstats(config, statfile, qcstats)
		barplot(config, qcstats['twsrefnorm'])
		removedupref(config)
	return statfile

def genomescan_run(config):
	counttablefile= config['genomescaninfo']['counttablefile'] if 'counttablefile' in config['genomescaninfo'] else os.path.join(config['resultdir'], 'counttable.txt.gz')
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
	makedirectory(outfile, config['dhmrinfo']['verbose'])
	hyperfile=outfile+'.hyper.bed'
	hypofile=outfile+'.hypo.bed'

	qcol = 6 if config['dhmrinfo']['method'] in ['ttest', 'chisq', 'gtest'] else 7
	cmd = "(printf '%%s\\n' chrom start end $(zcat %s | head -n 1 | cut -f 2-) | paste -s -d $'\\t'; zcat %s | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { getline } $%s<%s && $3>0 { n=split($1, a, \"[:-]\"); for (i=1; i<=n; i++) { printf a[i] OFS } for (j=2; j<=NF; j++) { printf $j ((j<NF)?OFS:ORS) } }' | sort -k 1,1 -k 2,2n -k 3,3n) > %s" % (
				testfile, testfile, qcol, config['dhmrinfo']['qthr'], hyperfile
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])
	cmd = "(printf '%%s\\n' chrom start end $(zcat %s | head -n 1 | cut -f 2-) | paste -s -d $'\\t'; zcat %s | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { getline } $%s<%s && $3<0 { n=split($1, a, \"[:-]\"); for (i=1; i<=n; i++) { printf a[i] OFS } for (j=2; j<=NF; j++) { printf $j ((j<NF)?OFS:ORS) } }' | sort -k 1,1 -k 2,2n -k 3,3n) > %s" % (
				testfile, testfile, qcol, config['dhmrinfo']['qthr'], hypofile
			)
	runcmdsh(cmd, echo=config['dhmrinfo']['verbose'])

	if config['dhmrinfo']['method'] in ['ttest', 'chisq', 'gtest']:
		runcmd("(head -n 1 %s; (bedtools merge -i %s -d %s -c 4,5,6,7,8 -o max,max,max,min,min -header | tail -n +2; bedtools merge -i %s -d %s -c 4,5,6,7,8 -o max,min,max,min,min -header | tail -n +2) | sort -k 8,8g ) | gzip -n > %s" % (
			hyperfile, hyperfile, config['dhmrinfo']['maxdistance'], hypofile, config['dhmrinfo']['maxdistance'], outfile
			), echo=config['dhmrinfo']['verbose'])
	else:
		runcmd("(head -n 1 %s; (bedtools merge -i %s -d %s -c 4,5,6,7,8,9 -o max,max,max,max,min,min -header | tail -n +2; bedtools merge -i %s -d %s -c 4,5,6,7,8,9 -o max,min,max,min,min,min -header | tail -n +2) | sort -k 9,9g ) | gzip -n > %s" % (
			hyperfile, hyperfile, config['dhmrinfo']['maxdistance'], hypofile, config['dhmrinfo']['maxdistance'], outfile
			), echo=config['dhmrinfo']['verbose'])
	runcmd('rm -f %s %s' % (hyperfile, hypofile), echo=config['dhmrinfo']['verbose'])

def dhmr_run(config, statfile, counttablefile):
	testfile = config['dhmrinfo']['testfile'] if 'testfile' in config['dhmrinfo'] else os.path.join(config['resultdir'], 'testfile.txt.gz')
	if not os.path.exists(testfile):
		makedirectory(testfile, config['dhmrinfo']['verbose'])
		if config['dhmrinfo']['method'] == 'ttest':
			ttest(config, statfile, counttablefile, testfile)
		elif config['dhmrinfo']['method'] == 'chisq':
			chisq(config, statfile, counttablefile, testfile)
		elif config['dhmrinfo']['method'] == 'gtest':
			gtest(config, statfile, counttablefile, testfile)
		elif config['dhmrinfo']['method'] == 'nbtest_sf':
			nbtest_sf(config, counttablefile, testfile)
		else:
			nbtest(config, statfile, counttablefile, testfile)
	dhmrfile = config['dhmrinfo']['dhmrfile'] if 'dhmrfile' in config['dhmrinfo'] else os.path.join(config['resultdir'], 'dhmr.txt.gz')
	if not os.path.exists(dhmrfile):
		mergedhmr(config, testfile, dhmrfile)

def inputfilter(config, counttablefile, testfile1, testfile2, inputfilterfile):
	keepregionfile=inputfilterfile + '.bed'
	cmd="(zcat %s | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { getline } $NF<%s && $3>0 { n=split($1, a, \"[:-]\"); for (i=1; i<=n; i++) { printf a[i] ((i<n)?OFS:ORS) } }'; zcat %s | awk -v FS='\\t' -v OFS='\\t' -e 'BEGIN { getline } $NF<%s && $3>0 { n=split($1, a, \"[:-]\"); for (i=1; i<=n; i++) { printf a[i] ((i<n)?OFS:ORS) } }') | sort -k 1,1 -k 2,2n -k 3,3n | uniq > %s" % (
			testfile1, config['inputinfo']['qthr'], testfile2, config['inputinfo']['qthr'], keepregionfile
			)
	runcmdsh(cmd, echo=config['inputinfo']['verbose'])
	runcmd("bedtools intersect -wa -a <(zcat %s) -b %s -header | gzip -n > %s" % (counttablefile, keepregionfile, inputfilterfile), echo=config['inputinfo']['verbose'])
	runcmd('rm -f %s' % (keepregionfile), echo=config['inputinfo']['verbose'])

import copy
def run(config):
	statfile=align_run(config)
	counttablefile=genomescan_run(config)
	if 'useinput' in config and config['useinput']:
		inputfilterfile=config['inputinfo']['inputfilterfile'] if 'inputfilterfile' in config['inputinfo'] else os.path.join(config['resultdir'], 'counttable_inputfilter.txt.gz')
		if not os.path.exists(inputfilterfile):
			testfile1=config['inputinfo']['testfile1'] if 'testfile1' in config['inputinfo'] else os.path.join(config['resultdir'], 'testfile_G1VsInput.txt.gz')
			dhmrfile1=config['inputinfo']['dhmrfile1'] if 'dhmrfile1' in config['inputinfo'] else os.path.join(config['resultdir'], 'testfile_G1VsInput.dhmr.gz')
			if not os.path.exists(testfile1):
				config_g1 = copy.deepcopy(config)
				config_g1['groupinfo']['group1']=config['groupinfo']['group1']
				config_g1['groupinfo']['group2']=config['inputinfo']['group1']
				config_g1['dhmrinfo']['method']=config['inputinfo']['method']
				config_g1['dhmrinfo']['testfile']=testfile1
				config_g1['dhmrinfo']['dhmrfile']=dhmrfile1
				dhmr_run(config_g1, statfile, counttablefile)
			testfile2=config['inputinfo']['testfile2'] if 'testfile2' in config['inputinfo'] else os.path.join(config['resultdir'], 'testfile_G2VsInput.txt.gz')
			dhmrfile2=config['inputinfo']['dhmrfile2'] if 'dhmrfile2' in config['inputinfo'] else os.path.join(config['resultdir'], 'testfile_G2VsInput.dhmr.gz')
			if not os.path.exists(testfile2):
				config_g2 = copy.deepcopy(config)
				config_g2['groupinfo']['group1']=config['groupinfo']['group2']
				config_g2['groupinfo']['group2']=config['inputinfo']['group2']
				config_g2['dhmrinfo']['method']=config['inputinfo']['method']
				config_g2['dhmrinfo']['testfile']=testfile2
				config_g2['dhmrinfo']['dhmrfile']=dhmrfile2
				dhmr_run(config_g2, statfile, counttablefile)
			inputfilter(config, config['genomescaninfo']['counttablefile'], testfile1, testfile2, inputfilterfile)
		counttablefile=inputfilterfile
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
		, formatter_class=argparse.RawDescriptionHelpFormatter
		, epilog='''
Example:
  cmsip -c cms.yaml

Date: 2020/06/10
Authors: Jin Li <lijin.abc@gmail.com>
'''
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
