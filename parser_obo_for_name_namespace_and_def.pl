#!/usr/bin/perl -w
use strict;
use POSIX;
use Data::Dumper;

## -- importing the script 'func.pl' depends on host OS, added on June 23, 2010 --
my $os = $^O;
if($os eq 'darwin'){ ## if mac --
    require '/Users/wchen/Dropbox/perl_scripts/func.pl';    
}elsif($os eq 'linux'){ ## -- if linux --
    require '/home/wchen/Dropbox/perl_scripts/func.pl';
}

## *********************************************
## ** version history **
## *********************************************
my $ver = '1.0';
my $last_modified = 'July 8, 2014';

## *********************************************
## ** GET opts **
## *********************************************
use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s");

if (!$opts{i} or !$opts{o} ){
    print "--------------------------------------------------------------------------------------------------
    \t\tversion : $ver by Weihua Chen; last modified : $last_modified
--------------------------------------------------------------------------------------------------
    USAGE: perl $0
        -i input go.obo file, all versions should do
        -o output file, go to name, namespace and defination 
--------------------------------------------------------------------------------------------------\n";
    exit;
}

## -- a liter version of obo-parser, available in file 'func.pl' --
my $refhash = &obo_parser_id2name_namespace_def( $opts{i} );
open OUT, ">$opts{o}";
print OUT join("\t", qw(goid name namespace def)), "\n";
while( my ( $acc, $arrref ) = each %{ $refhash } ){
    print OUT join("\t", $acc, @{ $arrref }), "\n";
}
close OUT;