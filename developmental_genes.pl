#!/usr/bin/perl -w
use strict;

## -- importing the script 'func.pl' depends on host OS, added on June 23, 2010 --
my $os = $^O;
my $go_children_script = "";
if($os eq 'darwin'){ ## if mac --
    require '/Users/wchen/Dropbox/perl_scripts/func.pl';
    $go_children_script = "/Users/wchen/Dropbox/perl_scripts/gene_ontology/go_children.pl"
}elsif($os eq 'linux'){ ## -- if linux --
    require '/home/wchen/Dropbox/perl_scripts/func.pl';
    $go_children_script = "/home/wchen/Dropbox/perl_scripts/gene_ontology/go_children.pl"
}

my $ver = '1.0';
my $last_modified = 'May 1, 2011';

## --------------------------------------------
## -- version history --
## --------------------------------------------

use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s","g:s","l:s","uid:i");

if (!$opts{i} or !$opts{o} or !$opts{g} ){
    print "----------------------------------------------------------------------
    \tver : $ver by Weihua Chen; last modified : $last_modified
    \tNOTE: developmental genes are those associated with two GO terms:
    \tGO:0007275\tGO:0030154
----------------------------------------------------------------------
    USAGE: perl $0
        -i input file contains a list of gene to go file [gene  GO:#####]
        -g input geneontology file, in obo 1.0 format
        -o output a list of developmental genes with GO annotation
            gene goID goTerm
      [optional]
        -l speicify to report genes in this file, default to report all genes
        -uid unique ids for current species, default is none
----------------------------------------------------------------------\n";
    exit;
}

## -- get genes from $opts{l} --
my %hValidGenes = ();
if( defined $opts{l} and -f $opts{l} ){
    %hValidGenes = &readfile($opts{l}, 0, 0);
}

## -- get goid to desc --
my $hrGOID2GOterm = &obo_parser_goID_to_name( $opts{g} );

## -- get a list of GO ids associated with 'development' --
my %hDevGOIDs = ();
my @aChildrenGOs = `perl $go_children_script -i 'GO:0007275' -i 'GO:0030154' -g $opts{g}`;
foreach my $line (@aChildrenGOs){
    if(defined $line and $line =~ /(GO:\d+)/){
        $hDevGOIDs{$1} ++;
    }
}
$hDevGOIDs{'GO:0030154'} ++;
$hDevGOIDs{'GO:0007275'} ++;

## -- check if input gene are developmental genes --
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
my %hDevGene2GO = ();
while(<IN>){
	next if( !/\S/ or /^#/); ## skip if current line is empty or contains annotation
    chomp;
    my ($gene, $go) = split(/\t/, $_);
    next if( !defined $go or !defined $gene or $go !~ /GO:\d+/);
    next if( !exists $hDevGOIDs{ $go }  or (defined $opts{l} and !exists $hValidGenes{$gene}) );
    my $goterm = exists $$hrGOID2GOterm{ $go } ?  $$hrGOID2GOterm{ $go } : "";
    $hDevGene2GO{ $gene }{ $go } = $goterm;
    
}
close IN;


open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
while( my ( $gene, $hrGO2Term ) = each %hDevGene2GO ){
    while(my ($goid, $goterm) = each %{$hrGO2Term}){
        print OUT $opts{uid}, "\t" if( defined $opts{uid} );
        print OUT join("\t", $gene, $goid, $goterm), "\n"
    }
}
close OUT;
