#!/usr/bin/perl -w
use strict;

# last modified : Feb 7, 2008
# By : chen
# HHUD

# find go term according to go id

use Getopt::Long;
#use Data::Dumper;

my %opts=();
GetOptions(\%opts,"i:s","o:s","g:s");

if (!$opts{i} or !$opts{o} or !$opts{g}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input file
        -g input gene_ontology.obo file
        -o out file
----------------------------------------------------------------------\n";
    exit;
}
my $start_time=time;

#### ontology parser
my %hNodesInfo = ();
my %hParents = ();
&obo_parser($opts{g},\%hNodesInfo,%hParents);

my %hGene2GOID = ();
#### process input native go file
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
$/="\n";
while(<IN>){
    chomp;
    next if (!/GO/);
    my ($geneID,@aGO) = split(/\t/,$_);
    
    foreach my $goid (@aGO){
        $hGene2GOID{$geneID}{$goid} = 1;
    }
}
close IN;
while (my ($geneID, $refHash)= each %hGene2GOID){
    foreach my $goID (keys %{$refHash}){
        if (exists $hNodesInfo{$goID} and !exists $hParents{$goID}){
            print OUT $geneID,"\t",$goID,"\t",$hNodesInfo{$goID}{name_space},"\t",$hNodesInfo{$goID}{name},"\n";
        }
    }
}
close OUT;

print "Time used : ",time-$start_time," seconds!!\n";

sub obo_parser{
    my ($infile,$hash,$parents)=@_;

    #### parsing obo_file
    my $backup=$/;
    $/="\n\[";
 
    my $header=1;
    open IN, $infile or die "Cannot open file: $infile!\n";
    while (<IN>){
        if ($header){
            $header=0;
            my ($ver)=/format-version:\s+([0-9.]+)/;
            print "\t===============================================================================\n\t\tParsing OBO file, current version is $ver .\n\t===============================================================================\n";
        }else{
            #### new [item]
            my $acc='';
            my $name_space='';
            my $name='';
            my $is_root=0;
            my $is_obsolete=0;
            
            my @aAlt_id=();
            my @aParents=();
            my @aRecommanded_ids_for_obsolete_id=();
            my @aContents=split(/\n/,$_);
            foreach my $line (@aContents){
                if ($line=~/^id: (\S+)/){
                    $acc=$1;
                    push @aAlt_id,$acc;
                }elsif($line=~/^alt_id: (\S+)/){
                    push @aAlt_id,$1;
                }elsif($line=~/^name: (.+)/){
                    $name=$1;
                }elsif($line=~/^namespace: (.+)/){
                    $name_space=$1;
                }elsif($line=~/^is_a:\s+(GO:\d+)/){
                    push @aParents,$1;
                }elsif($line=~/^relationship: part_of (GO:\d+)/){
                    push @aParents,$1;
                }elsif($line=~/is_obsolete: true/){
                    $is_obsolete=1;
                }elsif($line=~/^comment:/){
                    (@aRecommanded_ids_for_obsolete_id)=$line=~/(GO:\d+)/g;
                }
            }
            #### add new item in hash
            if ($is_obsolete){### if obsolete
                #foreach my $local_acc (@aAlt_id){
                #    $$obsoleteHash{$local_acc}=\@aRecommanded_ids_for_obsolete_id;
                #}
            }else{
                $is_root=1 if !(scalar @aParents);
                my %newHash=(acc=>$acc,parents=>\@aParents,name=>$name,name_space=>$name_space,is_root=>$is_root,is_relationship_build=>0);
                foreach my $local_acc (@aAlt_id){
                    $$hash{$local_acc}=\%newHash;
                }
                foreach my $parent (@aParents){
                    $$parents{$parent} = 1;
                }
            }
        }
    }
    close IN;
    
    #### building relationships for all nodes
    #foreach my $key (keys %{$hash}){
    #    if ((!$$hash{$key}{is_root}) and (!$$hash{$key}{is_relationship_build})){ # if not root and relationship not build
    #        for(my $i=0; $i<@{$$hash{$key}{parents}}; $i++){
    #            my $parent_go_id=$$hash{$key}{parents}[$i];
    #            $$hash{$key}{parents}[$i]=$$hash{$parent_go_id};
    #        }
    #        $$hash{$key}{is_relationship_build}=1;
    #    }
    #}
}