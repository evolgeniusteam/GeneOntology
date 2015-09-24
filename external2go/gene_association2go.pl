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
my $last_modified = 'Oct xx, 2014';

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
        -i input gene association file
        -o output gene2go file
      [note]
        format gene association file for a sigle organism for GO analyses
--------------------------------------------------------------------------------------------------\n";
    exit;
}

my %hash = ();
open IN, $opts{i} or die;
while(<IN>){
    chomp;
    next if( /^#/ or !/\S/ or /^!/);
    my ( undef, $id, undef, $go ) = split(/\s+/, $_);
    next if( $go !~ /GO/ );
    $hash{ $id }{ $go } ++;
}
close IN;

open OUT, ">$opts{o}" or die;
while( my ($id, $hrGos) = each %hash ){
    print OUT join("\t", $id, keys %{ $hrGos }), "\n";
}
close OUT;
