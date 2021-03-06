#!/usr/bin/perl -w
use strict;
use Getopt::Long;
#use Data::Dumper;

my %opts=();
GetOptions(\%opts,"i:s","o:s","l:s");

if (!$opts{i} or !$opts{o} or !$opts{l}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input ebi human uniprot gene associateion file
        -l input human gene list file
        -o out file
----------------------------------------------------------------------\n";
    exit;
    
}
my $start_time=time;

my %hIDtoGO = ();
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
while(<IN>){
    chomp;
    next if (/^!/); # skip annotation lines
    if(/\d+/){
        my ($db, $acc, $gene_id, $info, $GO_term, $goa, $evidence, $db_ref, $aspect, $annotation) = split(/\t/,$_);
        next if ($db ne 'UniProtKB');
        my ($alt_ids) = split(/:/,$annotation);
        $alt_ids =~ tr/ //d;
        my @aAltIDs = split(/,/,$alt_ids);
        my ($another_id) = split(/_/,$gene_id);
        push @aAltIDs, $acc, $gene_id, $another_id;
        foreach my $id (@aAltIDs){
            $hIDtoGO{$id}{$GO_term} = 1;
        }
    }
}
close IN;

open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
open LIST, $opts{l} or die "Cannot open file: $opts{l}!\n";
while(<LIST>){
    chomp;
    if(/^(\S+)/){
        my $geneID = $1;
        if (exists $hIDtoGO{$geneID}){
            print OUT $geneID,"\t", join("\t",keys %{$hIDtoGO{$geneID}}),"\n";
        }else{
            print $geneID,"\n";
        }
    }
}
close LIST;
close OUT;
print "Time used : ",time-$start_time," seconds!!\n";



