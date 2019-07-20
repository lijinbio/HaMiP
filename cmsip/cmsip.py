#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.0.1"

import os
import sys
import argparse
import yaml
import subprocess

def runcmd(cmd, log=subprocess.PIPE):
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

def bsmap_runcmd(fname, refenece, numthread, outfile):
	runcmd('mkdir -p ' + os.path.dirname(outfile))
	cmd = 'bsmap' + \
		' -a ' + fname + \
		' -d ' + refenece + \
		' -R -n 1 -r 0 ' + \
		' -p ' + str(numthread) + \
		' -o ' + outfile
	runcmd(cmd, log=open(outfile+".stdout", 'w+'))

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
			runcmd('samtools merge ' + outfile + ' ' + files)
		else:
			fname=os.path.join(config['datainfo']['fastqdir'], sampleinfo['filenames'][0])
			bsmap_runcmd(fname, config['datainfo'][reference], config['numthreads'], outfile)

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
	print('==>bsmap<==')
	# bsmap_ref(config, 'reference')
	# bsmap_ref(config, 'spikein')
	mpstat = {}
	mpstat['reference'] = bsmap_stat(config, 'reference')
	mpstat['spikein'] = bsmap_stat(config, 'spikein')
	return mpstat

def removeCommonReads_runcmd(infile1, infile2, outfile1, outfile2):
	bin=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'perl', 'removeCommonRead.pl')
	runcmd('mkdir -p ' + os.path.dirname(outfile1) + ' ' + os.path.dirname(outfile2))
	cmd = bin + ' ' + infile1 + ' ' + infile2 + \
		' >(' + 'samtools view -bS - ' + ' -o ' + outfile1 + ' 2>/dev/null)' \
		' >(' + 'samtools view -bS - ' + ' -o ' + outfile2 + ' 2>/dev/null)'
	cp=runcmd(cmd)
	return int(cp.stdout.strip())

def removeCommonReads(config):
	print('==>removeCommonReads<==')
	inbasedir=os.path.join(config['datainfo']['outdir'], 'bsmap')
	outbasedir=os.path.join(config['datainfo']['outdir'], 'removeCommonReads')
	comm={}
	for sampleinfo in config['sampleinfo']:
		refinfile=os.path.join(inbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkinfile=os.path.join(inbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		refoutfile=os.path.join(outbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkoutfile=os.path.join(outbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		comm[sampleinfo['sampleid']] = removeCommonReads_runcmd(refinfile, spkinfile, refoutfile, spkoutfile)
	return comm

def totalwigsums_n(f):
	cmd = "bedtools genomecov -ibam %s -bg | awk -e 'BEGIN { sum=0 } { sum += $4*($3-$2) } END { print sum }'" % f
	cp=subprocess.run(cmd, universal_newlines=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	if cp.returncode != 0:
		print('Error: %s failed.' % cmd, vars(cp), file=sys.stderr)
		sys.exit(-1)
	return int(cp.stdout.strip())

def totalwigsums(config):
	print('==>totalwigsums<==')
	inbasedir=os.path.join(config['datainfo']['outdir'], 'removeCommonReads')
	twss = {}
	for reference in ['reference', 'spikein']:
		tws = {}
		for sampleinfo in config['sampleinfo']:
			f=os.path.join(inbasedir, reference, sampleinfo['sampleid'] + '.bam')
			tws[sampleinfo['sampleid']] = totalwigsums_n(f)
		twss[reference] = tws
	return twss

import statistics
def estimateSizeFactors(tws):
	print('==>estimateSizeFactors<==')
	mws = statistics.median(tws.values())
	sizefactors = {id : mws / ws for id, ws in tws.items()}
	return sizefactors

def normalizetwsref(tws, sizefactors):
	print('==>normalizetwsref<==')
	twsn = {id: ws * sizefactors[id] for id, ws in tws.items()}
	return twsn

def removeDuplication():
	print('==>removeDuplication<==')

def normalizeWigsum():
	print('==>normalizeWigsum<==')

def DMR():
	print('==>DMR<==')

def run(config):
	mpstat = bsmap(config)
	print(mpstat)
	commcounts = removeCommonReads(config)
	print(commcounts)
	twss = totalwigsums(config)
	print(twss)
	sizefactors = estimateSizeFactors(twss['spikein'])
	print(sizefactors)
	twsrefnorm = normalizetwsref(twss['reference'], sizefactors)
	print(twsrefnorm)
	barplot(config, twsrefnorm)
	removeDuplication()
	normalizeWigsum()
	DMR()

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
