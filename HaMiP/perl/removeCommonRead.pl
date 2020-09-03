#!/usr/bin/env perl

use warnings;
use strict;
my ($f1, $f2, $fA, $fB)= @ARGV;
##print join("\t", $f1, $f2, $fA, $fB); exit;

my %read1 = ();
my %read2 = ();
my %comm = ();

open(READ1, "samtools view $f1 |") or die;
while(<READ1>){
	chomp;
	my @f=split(/\t/,$_);
	$read1{$f[0]}=1; 
}
close(READ1);


open(READ2, "samtools view $f2 |") or die;
while(<READ2>){
	chomp;
	my @f=split(/\t/,$_);
	if(defined $read1{$f[0]})
	{
		$comm{$f[0]}=1;
 	} else {
		$read2{$f[0]}=1;
	}
}
close(READ2);

open(READA, ">$fA") or die;
open(READ1, "samtools view -h $f1 |") or die;
while(<READ1>){
        chomp;
        my @f=split(/\t/,$_);
        if(defined $comm{$f[0]}){
	} else {
		print READA $_ , "\n";
	}
}
close(READ1);
close(READA);

open(READB, ">$fB") or die;
open(READ2, "samtools view -h $f2 |") or die;
while(<READ2>){
        chomp;
        my @f=split(/\t/,$_);
        if(defined $comm{$f[0]}){
        } else {
                print READB $_ , "\n";
        }
}
close(READ2);
close(READB);

my $size = keys %comm;
print $size, "\n";

