#!/usr/bin/perl -w
use strict;
use Getopt::Long;

## -- mac version
require '/Users/wchen/workspace/tools_perlscripts/func.pl';


my %opts=();
GetOptions(\%opts,"i:s","o:s", "r:s");

if (!$opts{i} or !$opts{o} or !$opts{r}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input gene list '/Users/wchen/workspace/databases/yeast/wang_duplicated_pair_and_kaks_and_age_and_expression/yeast_dup_versus_enssential_oct2009.tab'
	-r input gene to GO file '/Users/wchen/workspace/databases/yeast/SGD/literature_curation/go_slim_mapping.tab'
        -o out gene to GO
----------------------------------------------------------------------\n";
    exit;
}

## ===================================
## -- load input gene list --
my %hGene2Info = &readfile($opts{i}, 2 , 1);

## -- load gene to go list --
my %hGene2GO = ();
&gene2GO_parser_yeast($opts{r}, \%hGene2GO);

## -- print out --
open OUT, ">$opts{o}" or die;
while(my ($locus, $arrayref) = each %hGene2Info){
    my @aGO = ();
    if(exists $hGene2GO{$locus}){
	@aGO = sort keys %{$hGene2GO{$locus}};
    }
    print OUT join("\t", $locus, @{$arrayref}, @aGO), "\n";
}
close OUT;

#######################################
## -- sub functions --
sub gene2GO_parser_yeast{
    my ($infile, $return_hashref) = @_;
    open IN, $infile or die;
    while(<IN>){
	my ($locus, $gene, $serial, $cat, $name_space, $go, $type) = split(/\t/, $_);
	if($go =~ /GO/){
	    $$return_hashref{$locus}{$go} ++;
	    $$return_hashref{$gene}{$go} ++;	    
	}

    }
    close IN;
}