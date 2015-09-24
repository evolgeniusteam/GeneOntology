#!/usr/bin/perl -w
use strict;
use POSIX;
use Data::Dumper;
use String::Util qw(trim);
use Term::ANSIColor;


## -- importing the script 'go_functions.pl' depends on host OS, added on Sep 24, 2015 --
my $os = $^O;
if($os eq 'darwin'){ ## if mac --
    require '/Users/wchen/Dropbox/perl_scripts/gene_ontology/go_functions.pl';
}elsif($os eq 'linux'){ ## -- if linux --
    require '/home/wchen/Dropbox/perl_scripts/gene_ontology/go_functions.pl';
}

## *********************************************
## ** version history **
## *********************************************
my $ver = '1.0';
my $last_modified = 'Sep 24, 2015';

## *********************************************
## ** GET opts **
## *********************************************
use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s", "g:s","a:s");

if (!$opts{i} or !$opts{o}  or !$opts{g} or !$opts{a} ){
    print "--------------------------------------------------------------------------------------------------
    \t\tversion : $ver by Weihua Chen; last modified : $last_modified
--------------------------------------------------------------------------------------------------
    USAGE: perl $0
        -i input gene to go annotation file, it can be one of the following formats:
            gene    go ## a simple two-column file or
            annotation file downloaded from geneontology.com # at least 6 columns
             in which the 2nd column is gene name and the 6th column is the GO:####
        -g go.obo file
        -o output go to gene association results, it has two columns:
            go
            gene
        -a output go to name space to annotation file, it has three columns:
            go
            name_space
            description
--------------------------------------------------------------------------------------------------\n";
    exit;
}

## -- first, parse go.obo --
print STDERR "\tparsing go.obo file ... \n";
my %hGoObo = ();
my %hObsoleteTerms = ();
&obo_parser( $opts{g}, \%hGoObo, \%hObsoleteTerms);

## -- then, load gene to go file --
## -- $hash{ $go }{ $genes } ++
print STDERR "\tloading GO annotation file ... \n";
my $hrGo2Genes = &parse_geneontology_annotation_file( $opts{i} );

## -- get GO to gene file, a complete list;
## -- in this list, if a gene belongs to a GO, it should also belong to all the latter's parent GOs --
print STDERR "\tprocessing ...\n";
my %hGO2Genes = ();
while( my ($go, $hrGenes) = each %{ $hrGo2Genes } ){
    my $go_obj = $hGoObo{$go};
    next if( !defined $go_obj );

    ## -- check if root
    next if( $$go_obj{is_root} );
    &addMember2GO( $go_obj, $hrGenes, \%hGO2Genes );
}

print STDERR "\twrite to go -> gene association file ... \n";
open OUT, ">$opts{o}" or die;
while(my ( $go, $href ) = each %hGO2Genes){
    foreach my $gene (keys %{ $href }){
        print OUT join("\t", $go, $gene), "\n";
    }
}
close OUT;

print STDERR "\twrite to go -> name_space -> description file ... \n";
open OUT2, ">$opts{a}" or die;
while( my ($goid, $hashref) = each %hGoObo ){
    my $name_space = $$hashref{name_space};
    my $name = $$hashref{name};
    print OUT2 join("\t", $goid, $name_space, $name), "\n";
}
close OUT2;
print STDERR "\tall jobs done!!\n\n";

###################################
### sub functions --
## -- input : a GO node; (hash structure of obo_parser output { acc=..., parents = ...} )
##            a hashref to genes $hash{ $gene } +;
##            a hashref to hold final results
sub addMember2GO{
    my ( $go_obj, $hrGenes, $hrResults ) = @_;
    if( !$$go_obj{is_root} ){
        my $goid = $$go_obj{acc};
        my $go_aspect = $$go_obj{name_space};
        foreach my $gene ( keys %{ $hrGenes } ){
            $$hrResults{ $goid }{ $gene } = $go_aspect;
        }

        ## -- if hash parents
        foreach my $parent ( @{ $$go_obj{parents} } ){
            &addMember2GO( $parent, $hrGenes, $hrResults );
        }
    }
}


## -- INPUT: it can be one of the following formats:
#    gene    go ## a simple two-column file or
#    annotation file downloaded from geneontology.com # at least 6 columns
#     in which the 2nd column is gene name and the 6th column is the GO:####
## -- OUTPUT: $hash{ $go }{ $geneID } ++;
sub parse_geneontology_annotation_file{
    my ( $infile ) = @_;

    my %hash = ();
    open IN, $infile or die;
    while(<IN>){
        chomp;
        my @arr = split(/\t/, $_);

        my $gene;
        my $go;
        if ( scalar @arr == 2 ) {
            ( $gene, $go ) = @arr;
        } elsif ( scalar @arr >= 5 ) {
            ( undef, $gene, undef, undef, $go ) = @arr;
        }

        ## skip if gene or go IDs are not valid ...
        next if( !defined $gene or !defined $go or $go !~ /GO:\d+/ );
        $hash{ $go }{ $gene } ++;

    }
    close IN;
    return \%hash;
}
