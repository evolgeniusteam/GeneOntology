#!/usr/bin/perl -w
use strict;

use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input  file
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}
my %hGeneOntology=();
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while(<IN>){
    if (/^GR/){
        my @arr=split(/\t/,$_);
        #next if ($arr[12] ne 'taxon:39947');### japonica
        #print $arr[1],"\t",$arr[4],"\t",$arr[12],"\n";
        $hGeneOntology{$arr[1]}{$arr[4]}=1;
    }
}
close IN;
while(my ($acc,$ref)=each %hGeneOntology){
    print OUT $acc,"\t",join("\t", keys %{$ref}),"\n";
}
close OUT;
