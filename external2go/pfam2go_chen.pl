#!/usr/bin/perl -w
#Wei-Hua Chen
#2006-03-23

use strict;

use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s","m:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input pfam.parse file
        -o out file
	-d gene_association database, default PFAM, optional
        -m pfam id to acc mapper file, optonal  
----------------------------------------------------------------------\n";
    exit;
    
}

my $database='/database/GO/pfam2go';
if (defined $opts{d}){
	$database=$opts{d};
}

my $mapper='/database/Pfam_ls/id2acc.txt';
if (defined $opts{m}){
	$database=$opts{m};
}

my %hAcc2GO=();
open PFAM, "$database" or die "Cannot open database at $database!!\n";
while (<PFAM>){
    chomp;
    if (/^Pfam\:(\S+).*\s+(GO\:\d+)$/){
        if (exists $hAcc2GO{$1}){
            $hAcc2GO{$1}.=$2."\t";
        }else{
            $hAcc2GO{$1}=$2."\t";
        }
    }
}
close PFAM;


my %hID2GO=();
open MAPPER, "$mapper" or die "Cannot open database at $mapper!!\n";
while (<MAPPER>){
    chomp;
    if (my ($name,$acc)=/^(\S+)\t(PF\d+)/){
        if (exists $hAcc2GO{$acc}){
            $hID2GO{$name}=$hAcc2GO{$acc};
        }
    }
}
close MAPPER;

my %hPfam2GO=();
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
while(<IN>){
    chomp;
    next if (!/\d+/);
    if (my ($unigene_id,$domain)=/^(\S+)\s+(\S+)/){
        if (exists $hID2GO{$domain}){
            if (exists $hPfam2GO{$unigene_id}){
                $hPfam2GO{$unigene_id}.=$hID2GO{$domain};
            }else{
                $hPfam2GO{$unigene_id}=$hID2GO{$domain};
            }
        }
    }
}
close IN;

open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while (my ($k,$h)=each %hPfam2GO) {
	print OUT "$k\t$h\n";
}
close OUT;
