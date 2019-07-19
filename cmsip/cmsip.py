#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.0.1"

import argparse
import yaml
import os
import subprocess

def runcmd(cmd):
	print(cmd)
	print(subprocess.check_output(cmd, universal_newlines=True, shell=True))

def bsmap_runcmd(fname, refenece, numthread, outfile):
	runcmd('mkdir -p ' + os.path.dirname(outfile))
	cmd = 'bsmap' + \
		' -a ' + fname + \
		' -d ' + refenece + \
		' -R -n 1 -r 0 ' + \
		' -p ' + str(numthread) + \
		' -o ' + outfile
	runcmd(cmd)

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

def bsmap(config):
	print('==>bsmap<==')
	bsmap_ref(config, 'reference')
	bsmap_ref(config, 'spikein')

def removeCommonReads_runcmd(infile1, infile2, outfile1, outfile2):
	bin=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'perl', 'removeCommonRead.pl')
	runcmd('mkdir -p ' + os.path.dirname(outfile1) + ' ' + os.path.dirname(outfile2))
	# cmd = bin + ' ' + infile1 + ' ' + infile2 + \
	# 	' >(' + 'samtools view -bS - ' + ' -o ' + outfile1 + ' 2>/dev/null)' \
	# 	' >(' + 'samtools view -bS - ' + ' -o ' + outfile2 + ' 2>/dev/null)'
	# runcmd(cmd)
	sam1=os.path.splitext(outfile1)[0] + '.sam'
	sam2=os.path.splitext(outfile2)[0] + '.sam'
	cmd = bin + ' ' + infile1 + ' ' + infile2 + ' ' + sam1 + ' ' + sam2
	runcmd(cmd)
	runcmd('samtools view -bS ' + sam1 + ' ' + ' -o ' + outfile1 + ' 2>/dev/null')
	runcmd('samtools view -bS ' + sam2 + ' ' + ' -o ' + outfile2 + ' 2>/dev/null')
	runcmd('rm -f ' + sam1 + ' ' + sam2)

def removeCommonReads(config):
	print('==>removeCommonReads<==')
	inbasedir=os.path.join(config['datainfo']['outdir'], 'bsmap')
	outbasedir=os.path.join(config['datainfo']['outdir'], 'removeCommonReads')
	for sampleinfo in config['sampleinfo']:
		refinfile=os.path.join(inbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkinfile=os.path.join(inbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		refoutfile=os.path.join(outbasedir, 'reference', sampleinfo['sampleid'] + '.bam')
		spkoutfile=os.path.join(outbasedir, 'spikein', sampleinfo['sampleid'] + '.bam')
		removeCommonReads_runcmd(refinfile, spkinfile, refoutfile, spkoutfile)

def estimateSizeFactors():
	print('==>estimateSizeFactors<==')

def normalizeTotalWigsum():
	print('==>normalizeTotalWigsum<==')

def removeDuplication():
	print('==>removeDuplication<==')

def normalizeWigsum():
	print('==>normalizeWigsum<==')

def DMR():
	print('==>DMR<==')

def run(config):
	# bsmap(config)
	removeCommonReads(config)
	estimateSizeFactors()
	normalizeTotalWigsum()
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
