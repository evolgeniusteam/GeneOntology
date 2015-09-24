#!/usr/bin/perl -w
use strict;

## -- contains sub functions for GO analysis --


## INPUT: ogo file, hash_ref to hold results, hash_ref to hold obsolete GOs --
##        OUTPUT : 'results' is a hash that contains:
##                           $hash{ $go_acc  } =
                                # (
                                #     acc                   => $acc,
                                #     parents               => \@aParents, ## -- here parents contains [  {acc ... }, {}  ]
                                #     name                  => $name,
                                #     name_space            => $name_space,
                                #     is_root               => $is_root,
                                #     is_relationship_build => 0 ## will be 1 at the end of the run ...
                                # );
## depends on : none --
sub obo_parser{
    my ($infile, $hash, $obsoleteHash, $pureleaf_hashref, $return_altid_hashref)=@_;

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
    my ($goid, $hrNodesAndRelationship)=@_;
    my @aPath=();

    #print "goid is $goid\n";
    if (exists $$hrNodesAndRelationship{$goid}){
        my @aNewPath=($$hrNodesAndRelationship{$goid});
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

1;
