#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.0.1"

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
			print('Error: %s failed.' % cmd, vars(cp), file=sys.stderr)
			sys.exit(-1)
	except OSError as e:
		print("Execution failed: ", e, file=sys.stderr)
		sys.exit(-1)
	return cp

def bsmap_runcmd(fname, refenece, numthread, outfile, verbose=False):
	runcmd('mkdir -p ' + os.path.dirname(outfile), echo=verbose)
	cmd = 'bsmap' + \
		' -a ' + fname + \
		' -d ' + refenece + \
		' -R -n 1 -r 0 ' + \
		' -p ' + str(numthread) + \
		' -o ' + outfile
	runcmd(cmd, log=open(outfile+".stdout", 'w+'), echo=verbose)

def bsmap_ref(config, reference):
	outbasedir=os.path.join(config['datainfo']['outdir'], 'bsmap', reference)
	for sampleinfo in config['sampleinfo']:
		outfile=os.path.join(outbasedir, sampleinfo['sampleid'] + '.bam')
		if len(sampleinfo['filenames']) > 1:
			files = ''
			for f in sampleinfo['filenames']:
				bname=os.path.splitext(os.path.splitext(os.path.basename(f))[0])[0];
				fname=os.path.join(config['datainfo']['fastqdir'], f)
				singlefile=os.path.join(outbasedir, 'single', bname + '.bam')
				bsmap_runcmd(fname, config['datainfo'][reference], config['numthreads'], singlefile)
				files += ' ' + singlefile
			runcmd('samtools merge ' + outfile + ' ' + files, echo=config['verbose'])
		else:
			fname=os.path.join(config['datainfo']['fastqdir'], sampleinfo['filenames'][0])
			bsmap_runcmd(fname, config['datainfo'][reference], config['numthreads'], outfile, config['verbose'])

import re
def bsmap_stat_parse(infile):
	with open(infile) as f:
		dstr=f.read()
	totalreads=int(re.search('total reads: (\d+)', dstr).groups()[0])
	alignedreads=int(re.search('aligned reads: (\d+)', dstr).groups()[0])
	uniquereads=int(re.search('unique reads: (\d+)', dstr).groups()[0])
	return (totalreads, alignedreads, uniquereads)

def bsmap_stat(config, reference):
	basedir=os.path.join(config['datainfo']['outdir'], 'bsmap')
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
	if config['verbose']:
		print('==>bsmap<==')
	# bsmap_ref(config, 'reference')
	# bsmap_ref(config, 'spikein')
	mpstat = {}
	mpstat['reference'] = bsmap_stat(config, 'reference')
	mpstat['spikein'] = bsmap_stat(config, 'spikein')
	if config['verbose']:
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
	if config['verbose']:
		print('==>removeCommonReads<==')
	inbasedir=os.path.join(config['datainfo']['outdir'], 'bsmap')
	outbasedir=os.path.join(config['datainfo']['outdir'], 'removeCommonReads')
	comm={}
	for sampleinfo in config['sampleinfo']:
		refinfile=os.path.join(inbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkinfile=os.path.join(inbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		refoutfile=os.path.join(outbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkoutfile=os.path.join(outbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		comm[sampleinfo['sampleid']] = removeCommonReads_runcmd(refinfile, spkinfile, refoutfile, spkoutfile, config['verbose'])
	if config['verbose']:
		print(comm)
	return comm

def totalwigsums_n(f):
	cmd = "bedtools genomecov -ibam %s -bg | awk -e 'BEGIN { sum=0 } { sum += $4*($3-$2) } END { print sum }'" % f
	cp=subprocess.run(cmd, universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	if cp.returncode != 0:
		print('Error: %s failed.' % cmd, vars(cp), file=sys.stderr)
		sys.exit(-1)
	return int(cp.stdout.strip())

def totalwigsums(config):
	if config['verbose']:
		print('==>totalwigsums<==')
	inbasedir=os.path.join(config['datainfo']['outdir'], 'removeCommonReads')
	twss = {}
	for reference in ['reference', 'spikein']:
		tws = {}
		for sampleinfo in config['sampleinfo']:
			f=os.path.join(inbasedir, reference, sampleinfo['sampleid'] + '.bam')
			tws[sampleinfo['sampleid']] = totalwigsums_n(f)
		twss[reference] = tws
	if config['verbose']:
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

import matplotlib.pyplot as plt
def barplot(config, tws):
	outfile=os.path.join(config['datainfo']['outdir'], 'qcstats_twsn_barplot.pdf')
	plt.figure(figsize=(5, 5))
	ax=plt.axes()
	ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, loc: '{:,}'.format(int(x))))
	plt.bar(*zip(*tws.items()), width=0.5, color='blue')
	plt.title('Normalized total wigsum')
	runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['verbose'])
	plt.savefig(outfile, bbox_inches='tight')

def QCstats(config, qcstats):
	statfile=os.path.join(config['datainfo']['outdir'], 'qcstats.txt')
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

def removedupref(config):
	if config['verbose']:
		print('==>removedupref<==')
	indir=os.path.join(config['datainfo']['outdir'], 'removeCommonReads', 'reference')
	outdir=os.path.join(config['datainfo']['outdir'], 'removedupref')
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(indir, sampleinfo['sampleid'] + '.bam')
			outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bam')
			runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['verbose'])
			runcmd('samtools rmdup %s %s' % (infile, outfile), echo=config['verbose'])

import pandas as pd
def genomecov(config, statfile):
	if config['verbose']:
		print('==>genomecov<==')
	sizefactors=pd.read_csv(statfile, sep='\t', index_col='sample_id', usecols=['sample_id', 'sizefactors']).to_dict()['sizefactors']
	indir=os.path.join(config['datainfo']['outdir'], 'removedupref')
	outdir=os.path.join(config['datainfo']['outdir'], 'genomecov')
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(indir, sampleinfo['sampleid'] + '.bam')
			outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bedgraph')
			cmd = "bedtools genomecov -ibam %s -bg -scale %f > %s" % (infile, sizefactors[sampleinfo['sampleid']], outfile)
			runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['verbose'])
			runcmd(cmd, echo=config['verbose'])

def genomemeancov(config):
	if config['verbose']:
		print('==>genomemeancov<==')
	indir=os.path.join(config['datainfo']['outdir'], 'genomecov')
	outdir=os.path.join(config['datainfo']['outdir'], 'genomemeancov')
	for sampleinfo in config['sampleinfo']:
			infile=os.path.join(indir, sampleinfo['sampleid'] + '.bedgraph')
			outfile=os.path.join(outdir, sampleinfo['sampleid'] + '.bedgraph')
			cmd = "bedtools intersect -a %s -b %s -wo | awk -v OFS='\t' -e '{ print $1, $2, $3, $7*$8/%d }' | bedtools groupby -g 1,2,3 -c 4 -o sum > %s" % (config['datainfo']['windowfile'], infile, config['datainfo']['windowsize'], outfile)
			runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['verbose'])
			if config['verbose']:
				print(cmd)
			cp=subprocess.run(cmd, universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
			if cp.returncode != 0:
				print('Error: %s failed.' % cmd, vars(cp), file=sys.stderr)
				sys.exit(-1)

def meancovtable(config):
	if config['verbose']:
		print('==>meancovtable<==')
	indir=os.path.join(config['datainfo']['outdir'], 'genomemeancov')
	outfile=os.path.join(config['datainfo']['outdir'], 'meancovtable.txt.gz')
	sampleids=[sampleinfo['sampleid'] for sampleinfo in config['sampleinfo']]
	fs=[os.path.join(indir, id+'.bedgraph') for id in sampleids]
	cmd = "bedtools unionbedg -i %s -header -names %s | awk -v FS='\t' -v OFS='\t' -v separator='-' -v col=3 -e '{ for (i=1; i<col; i++) { printf $i ((i<NF)?separator:ORS) } for (i=col; i<=NF; i++) { printf $i ((i<NF)?OFS:ORS) } }' | gzip -n > %s" % (' '.join(fs), ' '.join(sampleids), outfile)
	runcmd('mkdir -p ' + os.path.dirname(outfile), echo=config['verbose'])
	if config['verbose']:
		print(cmd)
	cp=subprocess.run(cmd, universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	if cp.returncode != 0:
		print('Error: %s failed.' % cmd, vars(cp), file=sys.stderr)
		sys.exit(-1)

def DMR(config, cnttablefile):
	if config['verbose']:
		print('==>DMR<==')

def run(config):
	statfile = os.path.join(config['datainfo']['outdir'], 'qcstats.txt')
	if 'statfile' in config['datainfo']:
		statfile = config['datainfo']['statfile']
	if not os.path.exists(statfile):
		qcstats = {}
		qcstats['mpstat'] = bsmap(config)
		qcstats['comm'] = removeCommonReads(config)
		qcstats['twss'] = totalwigsums(config)
		qcstats['sizefactors'] = estimateSizeFactors(qcstats['twss']['spikein'], config['verbose'])
		qcstats['twsrefnorm'] = normalizetwsref(qcstats['twss']['reference'], qcstats['sizefactors'], config['verbose'])
		QCstats(config, qcstats)
		barplot(config, qcstats['twsrefnorm'])
		
	cnttablefile=os.path.join(config['datainfo']['outdir'], 'meancovtable.txt.gz')
	if 'cnttablefile' in config['datainfo']:
		cnttablefile = config['datainfo']['cnttablefile']
	if not os.path.exists(cnttablefile):
		if not os.path.exists(config['datainfo']['windowfile']):
			print('Error: the window file should be specified in the configuration file.'
					, 'Window file can be generated using bedtools makewindows.'
					, '  e.g. bedtools makewindows -g <(fetchChromSizes hg38) -w 100 > hg38_w100.bed'
					, sep='\n')
			sys.exit(-1)
		removedupref(config)
		genomecov(config, statfile)
		genomemeancov(config)
		meancovtable(config)
	DMR(config, cnttablefile)

def main():
	parser = argparse.ArgumentParser(
		description='CMS-IP sequencing analysis workflow'
		)
	parser.add_argument('-c', '--config'
		, type=argparse.FileType('r')
		, required=True
		, help='Configuration file in YAML format.'
		)
	args = parser.parse_args()
	config=yaml.load(args.config, Loader=yaml.FullLoader)
	run(config)

if __name__ == "__main__":
	main()
