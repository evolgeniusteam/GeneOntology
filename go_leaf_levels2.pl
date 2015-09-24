#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Data::Dumper;

my %opts=();
GetOptions(\%opts,"g:s","o:s","level:n", "namespace:s");

my $ver='1.0';
$ver = '2.0';
$ver = '3.0';

###########################################
## -- version history --
## -- version 2.0, Nov 4, 2009 --
## -- add support to 'alt_ids', which were missed in previous version --
## -- fix a bug in identifying 'pure-leaf' nodes, which is very serious; however, this bug didn't change many of the results, which is strange ...
## -- now the 'DistanceFromLeaf' stands for the shortest path of current node to any leaves --
###########################################

if (!$opts{g} or !$opts{o}){
    print "========================================================================
                    Current version is $ver
========================================================================
    USAGE: perl $0
        -g input obo_go file
	-o output GO_terms (GO:#######) at levels up
      [optional]
        -level, 0 to 20, default is 2
	*** note: suppose you have a leaf GOID, it goes a long way to root 'the ontology'
	    then, the path will be numbered as the following:
	    0       1               2           3         4          5
	    leaf -> internal note1 -> Inode2 -> Inode3 -> Inode 4 -> root
	    if you choose 2, for example, GOids at position 0, 1 will be extracted
	*** NOTE: this script only works well on OBO format Ver1.0 NOT 1.2 ****
	-namespace, molecular_function|biological_process, default is all
========================================================================\n";
    exit;
}

## ==================================================
## -- global variables --
my $level=defined $opts{level} ? $opts{level} : 2;
$level = 2 if($level < 0);

#### ontology parser
my %hNodesAndRelationship=();
my %hObsoleteTerms=();
my %hPureLeafNodes = ();
my %hAltIDs = (); ## -- hash{master_GO_id}{$alterids} ++;
&obo_parser($opts{g},\%hNodesAndRelationship,\%hObsoleteTerms,\%hPureLeafNodes, \%hAltIDs);

## -- print Dumper %hPureLeafNodes,"\n";
#### find internal node close to the leaf node
my %hResultsTemp = ();
while(my ($key, $name_space) = each %hPureLeafNodes){
	next if( defined $opts{namespace} and $opts{namespace} ne $name_space );
    my $allpath_array_ref = &paths_to_root($key);
    foreach my $current_path (@{$allpath_array_ref}){
		my @aPath = reverse @{$current_path};
		for(my $i = 0; $i < $level; $i ++){
			if(defined $aPath[$i]){
				my $current_node = $aPath[$i]{acc};
				my @arr = ($i, $key, $name_space);
				push @{$hResultsTemp{$current_node}}, \@arr;
			}
		}
    }
}

## -- searching for the shortest path from current node to any leaf nodes --
my %hResults = ();
while(my ($key, $refarray) = each %hResultsTemp){
    my @arr = sort {$a->[0] <=> $b->[0]} @{$refarray};
    $hResults{$key} = $arr[0];
}

## --  output --
open OUT, ">$opts{o}" or die;
print OUT join("\t", qw(GO LevelFromLeaf Leaf NameSpace)),"\n";
while(my ($key, $h) = each %hResults){
    ## -- master ID --
    print OUT join("\t", $key, @{$h}),"\n";
    
    ## -- alternative IDs --
    if(exists $hAltIDs{$key}){
		foreach my $altID (keys %{$hAltIDs{$key}}){
			print OUT join("\t", $altID, @{$h}), "\n";
		}
    }
}
close OUT;


sub obo_parser{
    my ($infile,$hash,$obsoleteHash,$pureleaf_hashref, $return_altid_hashref)=@_;

    #### parsing obo_file
    ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ## -- note: this part is not the same to its orginal version in go_chenwh_ver#.pl
    
    my $backup=$/;
    $/="\n\[";
 
    my $header=1;
    open IN, $infile or die "Cannot open file: $infile!\n";
    while (<IN>){
        if ($header){
            $header=0;
            my ($ver)=/format-version:\s+([0-9.]+)/;
            print STDERR "\t===============================================================================\n\t\tParsing OBO file, format ver $ver .\n\t===============================================================================\n";
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
                }elsif($line=~/^alt_id: (\S+)/){
                    $$return_altid_hashref{$acc}{$1} ++;
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
		$$hash{$acc} = \%newHash;
            }
        }
    }
    close IN;
        
    my %hParentNodes = (); ## -- nodes that are parent to some other node --
    #### building relationships for all nodes
    foreach my $key (keys %{$hash}){
        if ((!$$hash{$key}{is_root}) and (!$$hash{$key}{is_relationship_build})){ # if not root and relationship not build
            for(my $i=0; $i<@{$$hash{$key}{parents}}; $i++){
                my $parent_go_id=$$hash{$key}{parents}[$i];
                $$hash{$key}{parents}[$i]=$$hash{$parent_go_id};
		$hParentNodes{$parent_go_id} ++; ## -- BUG here in version 1.0 --
            }
            $$hash{$key}{is_relationship_build}=1;
        }
    }
    
    foreach my $key (keys %{$hash}){
	if((!$$hash{$key}{is_root}) and !exists $hParentNodes{$key}){
	    $$pureleaf_hashref{$key} = $$hash{$key}{name_space} if($$hash{$key}{name_space} ne 'cellular_component'); ## -- cellular component are removed --
	}
    }
    
    print STDERR "\t\tthere are in total ", scalar keys  %{$hash}, " nodes in this go file ... \n";
    print STDERR "\t\tin which ", scalar keys %{$pureleaf_hashref}, " are leaf nodes ... \n\n";
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
