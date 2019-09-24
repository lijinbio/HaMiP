#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -v
../cmsip/cmsip.py -c cms.yaml -D genomescaninfo.readscount=True -D genomescaninfo.fragsize=200 -D dhmrinfo.qthr=0.05 -D resultdir=out_result groupinfo.group1=K groupinfo.group2=W
