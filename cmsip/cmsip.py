#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.0.1"

import argparse
import yaml
import os
import subprocess

def runcmd(cmd):
	print(subprocess.check_output(cmd, universal_newlines=True, shell=True))

def bsmap_runcmd(fname, refenece, numthread, outfile):
	runcmd('mkdir -p ' + os.path.dirname(outfile))
	cmd = 'bsmap' + \
		' -a ' + fname + \
		' -d ' + refenece + \
		' -n 1 -r 0 ' + \
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

def removeCommonReads():
	print('==>removeCommonReads<==')

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
	bsmap(config)
	removeCommonReads()
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
