#!/usr/bin/perl -w

#a siple ko2go mapping
#the input file should be like: 'GeneID koid'
#the output file is to be like: 'geneID GO:#####    GO:#####    GO:#####'
#Weihua Chen
#chenwh@genomics.org.cn
#2006-8-11

use strict;
use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input pfam.parse file
        -o out file
	-d gene_association database, default ko2go, optional
----------------------------------------------------------------------\n";
    exit;
    
}

#load database
my $database=defined $opts{d} ? $opts{d} : '/database/GO/ko2go';
my %hKO2GO=();
open KOGO, "$database" or die "Cannot open database at $database!!\n";
while (<KOGO>){
    chomp;
    next if (/^#/);
    if (my ($ko,$desc)=/^(\S+)\s+(.*)/){
        my $go='';
        while ($desc=~/(\d+)/g){
            $go.="\tGO:".$1;
        }
        if ($go){
            $hKO2GO{$ko}=$go;
        }
    }
}
close KOGO;

my $count=0;
my $ko_annotated_item=0;
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while(<IN>){
    chomp;
    if (my ($gene,$ko)=/^(\S+)\t(\S+)/){
        $ko_annotated_item++;
        if (exists $hKO2GO{$ko}){
            print OUT $gene,$hKO2GO{$ko},"\n";
            $count++;
        }
    }
}
close IN;
close OUT;
print "===========================================================================
	\t# of KO annotated items : $ko_annotated_item
        \t# of items mapped to go : $count
===========================================================================\n";
