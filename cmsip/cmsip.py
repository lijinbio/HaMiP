#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

__version__ = "0.0.0.1"

import argparse
import yaml

def bsmap(config):
	print('==>bsmap<==')

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
