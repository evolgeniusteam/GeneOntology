#!/usr/bin/perl -w
use strict;

my $ver = '1.0';
my $last_modified = 'May 3, 2011';

## -- importing the script 'func.pl' depends on host OS, added on June 23, 2010 --
my $os = $^O;
if($os eq 'darwin'){ ## if mac --
    require '/Users/wchen/Dropbox/perl_scripts/gene_ontology/go_lib.pl';    
}elsif($os eq 'linux'){ ## -- if linux --
    require '/home/wchen/Dropbox/perl_scripts/gene_ontology/go_lib.pl';
}

## --------------------------------------------
## -- version history --
## --------------------------------------------

use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s");

if (!$opts{i} or !$opts{o}){
    print "----------------------------------------------------------------------
    \t\t\tversion : $ver by Weihua Chen; last modified : $last_modified
----------------------------------------------------------------------
    USAGE: perl $0
        -i input go.obo file
        -o out go term to aspect file  
----------------------------------------------------------------------\n";
    exit;
}

#### ontology parser
my %hNodesAndRelationship=();
my %hObsoleteTerms=();
my %hPureLeafNodes = ();
my %hAltIDs = (); ## -- hash{master_GO_id}{$alterids} ++;
&obo_parser($opts{i},\%hNodesAndRelationship,\%hObsoleteTerms,\%hPureLeafNodes, \%hAltIDs);

open OUT, ">$opts{o}" or die;
print OUT join("\t", qw(GOTerm Name Aspect isLeaf)),"\n";
while( my ($k, $v) = each %hNodesAndRelationship ){
	next if( $k !~ /^GO/ or $$v{'is_root'} > 0 );
	print OUT join("\t", $k, $$v{'name'}, $$v{'name_space'}, exists $hPureLeafNodes{ $k } ? "Y" : "N" ), "\n";
}