#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my %opts=();
GetOptions(\%opts,"g:s","o:s","i:s@","r:s","debug","res:s");


my $ver='4.0';

if (!$opts{g} or !$opts{i} or !$opts{o} or !$opts{r}){
    print "========================================================================
                    Current version is $ver
=========================================================================
    USAGE: perl $0
        -i input native go file, GENEID GO:######, multiple, no more than 8
        -g input obo_go file
        -r input reference GO list
        -o print gene ontology result as text file
      [optional]
        -res input gene list file, if true, only print genes listed in this file
========================================================================\n";
    exit;
}
my %hRestrictedGeneList = ();
if(exists $opts{res}){
    open RES, $opts{res} or die;
    while(<RES>){
        chomp;
        if(/(\S+)/){
            $hRestrictedGeneList{$1} ++;
        }
    }
    close RES;   
}


my %hRefGO = ();
open REFGO, $opts{r} or die;
while(<REFGO>){
    chomp;
    if(/(GO:\d+)/){
        $hRefGO{$1} ++;
    }
}
close REFGO;

my %hGO2Genes=();
my %hResults = ();

my @nm_go_file=@{$opts{i}};
my $fileNumCount=scalar @nm_go_file;

my @maxNumOfUnigenes=();
my $maxNumGenes=0;

for (my $i=0; $i<@nm_go_file; $i++){
    my %hUniqGeneIDs = ();
    open IN, $nm_go_file[$i] or die "Cannot open file: $nm_go_file[$i]!\n";
    while(<IN>){
        next if (!/GO\:\d+/);
        chomp;
        my ($gene,@aGOIDs)=split(/\s+/,$_);
        $hUniqGeneIDs{$gene} = 1;
        
        #print Dumper (@aGOIDs),"\n";
        foreach my $goid (@aGOIDs){
            $hGO2Genes{$goid}{$i}{$gene}=1; ### add file_serial, starting at 0
        }
    }
    close IN;
    my $geneCount = scalar keys %hUniqGeneIDs;
    $maxNumGenes=$geneCount if ($geneCount>$maxNumGenes);
    push @maxNumOfUnigenes,$geneCount;
}

#### ontology parser
my %hNodesAndRelationship=();
my %hObsoleteTerms=();
&obo_parser($opts{g},\%hNodesAndRelationship,\%hObsoleteTerms);

## searching for putative obsolete go items ..
foreach my $goid (keys %hGO2Genes){
    if (exists $hObsoleteTerms{$goid}){
        my $refhash=$hGO2Genes{$goid};
        foreach my $alt_goid (@{$hObsoleteTerms{$goid}}){
            $hGO2Genes{$alt_goid}=$refhash;
        }
        delete $hGO2Genes{$goid};
    }
}

#### GO search by term
while (my ($goid,$refGeneHash)=each %hGO2Genes){
    print $goid,"\n" if(exists $opts{debug});
    my $allPath=&paths_to_root($goid);
    #print "now",Dumper(@{$allPath});
    my $i = 0;
    foreach my $local_path (@{$allPath}){
        $i ++;
        print "path # $i:\n" if(exists $opts{debug});
        my $j = 0;
        foreach my $level (@{$local_path}){
            $j ++;
            my $acc = $$level{acc};
            #print $acc,"\n";
            #print "\tnode $j:\t", $acc,"\n";
            if(exists $hRefGO{$acc}){
                print "exists : ", $acc,"\n" if(exists $opts{debug});
                while(my ($file_serial, $refhash_genes) = each %{$refGeneHash}){
                    foreach my $gene (keys %{$refhash_genes}){
                        $hResults{$acc}{genes}{$file_serial}{$gene} ++;
                    }
                }
                $hResults{$acc}{name} = $hNodesAndRelationship{$acc}{name};
            }
        }
    }
}

print Dumper(%hResults),"\n========================\n" if(exists $opts{debug});

open OUT, ">$opts{o}"  or die;
foreach my $GO (sort keys %hResults){
    my $name = $hResults{$GO}{name};
    print OUT join("\t", $GO, $name),"\n";
    foreach my $infile (sort{$a<=>$b} keys %{$hResults{$GO}{genes}}){
        my @aGenes = keys %{$hResults{$GO}{genes}{$infile}};
        if(exists $opts{res}){
            my @aRestrictedGeneList = ();
            foreach my $gene(@aGenes){
                push @aRestrictedGeneList, $gene if(exists $hRestrictedGeneList{$gene});
            }
            @aGenes = @aRestrictedGeneList;
        }
        print OUT "\t", join("\t", $nm_go_file[$infile], $maxNumOfUnigenes[$infile], scalar @aGenes, join(",", sort @aGenes)),"\n";
    }
}
close OUT;

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
            }
            $$hash{$key}{is_relationship_build}=1;
        }
    }
}

sub paths_to_root{
    my ($goid)=@_;
    my @aPath=();
    
    #print "goid is $goid\n";
    if (exists $hNodesAndRelationship{$goid}){
        my @aNewPath=($hNodesAndRelationship{$goid});
        push @aPath,\@aNewPath;
        my $has_parents=1;
        while ($has_parents){
            $has_parents=0;
            #### iterate all path
            foreach my $path (@aPath){
                #### check the top node for every path
                if (scalar @{$$path[0]{parents}}){
                    $has_parents++;
                    for (my $i=1; $i<@{$$path[0]{parents}}; $i++){
                        #### duplicate current path
                        my $newPath=&duplicate_path($path);
                        #### add new top node to deplicated path
                        unshift  @{$newPath},$$path[0]{parents}[$i];
                        push @aPath, $newPath;  
                    }
                    #### update current path
                    unshift  @{$path},$$path[0]{parents}[0];
                }else{#### no parents means current node top node
                    #### do nothing
                }
            }
        }
    }
    return \@aPath;
}

sub duplicate_path{
    my ($path)=@_;
    my @aNewPath=();
    foreach my $arrayElement (@{$path}){
        push @aNewPath,$arrayElement;
    }
    return \@aNewPath;
}

