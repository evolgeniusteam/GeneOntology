#!/usr/bin/perl -w
use strict;

use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input file
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}

my %hResults=();
my $cat='';
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
print OUT "GO level 1	GO level 2	#	%","\n";
while(<IN>){
    if (/<(.*)>/){
        $cat=$1;
    }elsif(/^  \%(.*)/){
        my ($term,$count,$percent)=split(/\t/,$1);
        chop $percent;
        $hResults{$cat}{$term}="$count\t$percent";
    }
}
close IN;
my @aCat=sort keys %hResults;
foreach my $cat (@aCat){
    print OUT $cat;
    my @aTerm=sort keys %{$hResults{$cat}};
    foreach my $term (@aTerm){
        print OUT "\t",$term,"\t",$hResults{$cat}{$term},"\n";
    }
}

close OUT;
#######
### Takes output of go_chenwh.pl as input,
###