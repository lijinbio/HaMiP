#!/usr/bin/env bash

function cmd {
local outfile=$1
local url=$2
[ -f "$outfile" ] || wget -q -O "$outfile" "$url"
}

mkdir -p fasta
cmd fasta/hg38.chr4.fa https://ndownloader.figshare.com/files/23080391?private_link=edab08f11290a321c507
cmd fasta/mm10.chr3.fa https://ndownloader.figshare.com/files/23135063?private_link=c9ecc9a289ae0032134b
