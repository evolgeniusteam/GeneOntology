#!/usr/bin/perl -w
use strict;
use Getopt::Long;


my $ver = '1.0';

use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s@", "g:s");

if (!$opts{i} or !$opts{g}){
    print "----------------------------------------------------------------------
    \t\t\tversion : $ver
----------------------------------------------------------------------
    USAGE: perl $0
        -i input GOID (GO:######), could be multiple, such as -i GO:#### -i GO:#####
        -g input go.obo file
      [note]
        this script will print the input GOID and all its children to stdout, one per line
----------------------------------------------------------------------\n";
    exit;
    
}

#### ontology parser
my %hNodesAndRelationship           = ();
my %hObsoleteTerms                  = ();
&obo_parser($opts{g},\%hNodesAndRelationship,\%hObsoleteTerms);


#use Data::Dumper;
#print Dumper %hNodesAndRelationship;
#exit;

my %hChildren                   = ();
foreach my $goid (@{$opts{i}}){
    &find_children($goid, $hNodesAndRelationship{$goid}, \%hChildren);
    
    $hChildren{$goid} ++;
}
print join("\n", sort keys %hChildren), "\n";

## --
sub find_children{
    my ($goid, $hrCurrentNode, $hrChildren)=@_;
    return if (!defined $$hrCurrentNode{children});
    
    foreach my $child_node (@{$$hrCurrentNode{children}}){
        my $current_goid            = $$child_node{acc};
        $$hrChildren{$current_goid} ++;
        
        &find_children($current_goid, $child_node, $hrChildren);
    }
}

sub obo_parser{
    my ($infile,$hash,$obsoleteHash)=@_;

    #### parsing obo_file
    my $backup=$/;
    $/="\n\[";
 
    my $header=1;
    open IN, $infile or die "Cannot open file: $infile!\n";
    while (<IN>){
        if ($header){
            $header=0;
            my ($ver)=/format-version:\s+([0-9.]+)/;
            print STDERR "\t===============================================================================\n\t\tParsing OBO file, current version is $ver .\n\t===============================================================================\n";
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
                foreach my $local_acc (@aAlt_id){
                    $$obsoleteHash{$local_acc}=\@aRecommanded_ids_for_obsolete_id;
                }
            }else{
                $is_root=1 if !(scalar @aParents);
                my %newHash=(acc=>$acc,parents=>\@aParents,name=>$name,name_space=>$name_space,is_root=>$is_root,is_relationship_build=>0);
                foreach my $local_acc (@aAlt_id){
                    $$hash{$local_acc}=\%newHash;
                }
            }
        }
    }
    close IN;
    
    #### building relationships for all nodes
    foreach my $key (keys %{$hash}){
        if ((!$$hash{$key}{is_root}) and (!$$hash{$key}{is_relationship_build})){ # if not root and relationship not build
            for(my $i=0; $i<@{$$hash{$key}{parents}}; $i++){
                my $parent_go_id=$$hash{$key}{parents}[$i];
                $$hash{$key}{parents}[$i]=$$hash{$parent_go_id};
                
                push @{$$hash{$parent_go_id}{children}}, $$hash{$key};
            }
            $$hash{$key}{is_relationship_build}=1;
        }
    }
}
