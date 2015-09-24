#!/usr/bin/perl -w

######## functions related to gene ontology analysis ####
use strict;


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
sub obo_parser {
    my ( $infile, $hash, $obsoleteHash ) = @_;

    #### parsing obo_file
    my $backup = $/;
    $/ = "\n\[";

    my $header = 1;
    open IN, $infile or die "Cannot open file: $infile!\n";
    while (<IN>) {
        if ($header) {
            $header = 0;
            my ($ver) = $_ =~ /format-version:\s+([0-9.]+)/;
            print STDERR "\t===============================================================================\n";
            print STDERR "\t\tParsing OBO file, current version is $ver .\n";
            print STDERR "\t===============================================================================\n";
        }
        else {
            #### new [item]
            my $acc         = '';
            my $name_space  = '';
            my $name        = '';
            my $is_root     = 0;
            my $is_obsolete = 0;

            my @aAlt_id                          = ();
            my @aParents                         = ();
            my @aRecommanded_ids_for_obsolete_id = ();
            my @aContents                        = split( /\n/, $_ );
            foreach my $line (@aContents) {
                if ( $line =~ /^id: (\S+)/ ) {
                    $acc = $1;
                    push @aAlt_id, $acc;
                }
                elsif ( $line =~ /^alt_id: (\S+)/ ) {
                    push @aAlt_id, $1;
                }
                elsif ( $line =~ /^name: (.+)/ ) {
                    $name = $1;
                }
                elsif ( $line =~ /^namespace: (.+)/ ) {
                    $name_space = $1;
                }
                elsif ( $line =~ /^is_a:\s+(GO:\d+)/ ) {
                    push @aParents, $1;
                }
                elsif ( $line =~ /^relationship: part_of (GO:\d+)/ ) {
                    push @aParents, $1;
                }
                elsif ( $line =~ /is_obsolete: true/ ) {
                    $is_obsolete = 1;
                }
                elsif ( $line =~ /^comment:/ ) {
                    (@aRecommanded_ids_for_obsolete_id) = $line =~ /(GO:\d+)/g;
                }
            }
            #### add new item in hash
            if ($is_obsolete) {    ### if obsolete
                foreach my $local_acc (@aAlt_id) {
                    $$obsoleteHash{$local_acc} =
                      \@aRecommanded_ids_for_obsolete_id;
                }
            }
            else {
                $is_root = 1 if !( scalar @aParents );
                my %newHash = (
                    acc                   => $acc,
                    parents               => \@aParents,
                    name                  => $name,
                    name_space            => $name_space,
                    is_root               => $is_root,
                    is_relationship_build => 0
                );
                foreach my $local_acc (@aAlt_id) {
                    $$hash{$local_acc} = \%newHash;
                }
            }
        }
    }
    close IN;

    $/ = $backup; 

    #### building relationships for all nodes
    foreach my $key ( keys %{$hash} ) {
        if (    ( !$$hash{$key}{is_root} )
            and ( !$$hash{$key}{is_relationship_build} ) )
        {    # if not root and relationship not build
            for ( my $i = 0 ; $i < @{ $$hash{$key}{parents} } ; $i++ ) {
                my $parent_go_id = $$hash{$key}{parents}[$i];
                $$hash{$key}{parents}[$i] = $$hash{$parent_go_id};
            }
            $$hash{$key}{is_relationship_build} = 1;
        }
    }
}

1;
