#!/usr/bin/env bash

cat <<EOF > config.yaml
sampleinfo:
  - sampleid: T1
    group: T
    reference:
      - $PWD/bamfile/T1.bam
    spikein:
      - $PWD/bamfile/T1_spk.bam
  - sampleid: T2
    group: T
    reference:
      - $PWD/bamfile/T2.bam
    spikein:
      - $PWD/bamfile/T2_spk.bam
  - sampleid: W1
    group: W
    reference:
      - $PWD/bamfile/W1.bam
    spikein:
      - $PWD/bamfile/W1_spk.bam
  - sampleid: W2
    group: W
    reference:
      - $PWD/bamfile/W2.bam
    spikein:
      - $PWD/bamfile/W2_spk.bam
  - sampleid: I1
    group: I
    reference:
      - $PWD/bamfile/I1.bam
    spikein:
      - $PWD/bamfile/I1_spk.bam
  - sampleid: I2
    group: I
    reference:
      - $PWD/bamfile/I2.bam
    spikein:
      - $PWD/bamfile/I2_spk.bam
groupinfo:
  group1: T
  group2: W
resultdir: $PWD/outdir
aligninfo:
  reference: $PWD/fasta/hg38.chr4.fa
  usespikein: True
  inputbam: True
  spikein: $PWD/fasta/mm10.chr4.fa
  statfile: $PWD/outdir/qcstats.txt
  barplotinfo:
    outfile: $PWD/outdir/qcstats_twsn_barplot.pdf
    height: 5
    width: 5
  numthreads: 4
  verbose: True
genomescaninfo:
  readextension: True
  fragsize: 100
  windowfile: $PWD/outdir/hg38_w100.bed
  referencename: hg38
  windowsize: 100
  readscount: True
  counttablefile: $PWD/outdir/extTrue_cntTrue_w100.txt.gz
  verbose: True
dhmrinfo:
  method: gtest
  meandepth: 1
  testfile: $PWD/outdir/gtest_extTrue_cntTrue_w100.txt.gz
  qthr: 0.05
  maxdistance: 0
  dhmrfile: $PWD/outdir/gtest_extTrue_cntTrue_w100.dhmr.gz
  numthreads: 4
  nsplit: 1000
  verbose: True
  keepNA: True
useinput: True
inputinfo:
  group1: I
  group2: I
  method: gtest
  qthr: 0.05
  testfile1: $PWD/outdir/gtest_extTrue_cntTrue_w100_TvsI.txt.gz
  dhmrfile1: $PWD/outdir/gtest_extTrue_cntTrue_w100_TvsI.dhmr.gz
  testfile2: $PWD/outdir/gtest_extTrue_cntTrue_w100_WvsI.txt.gz
  dhmrfile2: $PWD/outdir/gtest_extTrue_cntTrue_w100_WvsI.dhmr.gz
  inputfilterfile: $PWD/outdir/extTrue_cntTrue_w100_inputfilter.txt.gz
  verbose: True
EOF
