#!/usr/bin/perl -w
use strict;

my $ver = '1.0';
my $last_modified = 'Feb 19, 2011';

## --------------------------------------------
## -- version history --
## --------------------------------------------

use Getopt::Long;
my %opts=();
GetOptions(\%opts,"i:s","o:s", "j:s", "t:s");

if (!$opts{i} or !$opts{o}  or !$opts{j}  or !$opts{t}){
    print "----------------------------------------------------------------------
    \t\tversion : $ver by Weihua Chen; last modified : $last_modified
----------------------------------------------------------------------
    USAGE: perl $0
        -i input _se.txt file
        -j input go to gene file (_gce.txt)
        -t input type, all|under|over, default = all
	-o output file
----------------------------------------------------------------------\n";
    exit;
}

## ==== load go to gene file ===
my %hGO2genes = (); ## $hash{$go}{$type}{$gene} ++
open GOGENE, $opts{j} or die;
while(<GOGENE>){
    chomp;
    my ( $go, $cat, $gene, $type ) = split(/\t/, $_);
    if( $go =~ /GO:(\d+)/ ){
        my $goserial = $1;
        $hGO2genes{ $goserial }{ 'all' }{ $gene } ++;
        $hGO2genes{ $goserial }{ $type }{ $gene } ++;
    }
}
close GOGENE;


my $title = 1;
open OUT,">$opts{o}" or die "Cannot create file: $opts{o}!\n";
open IN, $opts{i} or die "Cannot open file: $opts{i}!\n";
while(<IN>){
    chomp;
    if( $title ){
	$title = 0;
        print OUT $_, "\t", "genes", "\n";
    } else {
        my ($goserial, @arr)  = split(/\t/, $_);
        my $genes = "";
        if( exists $hGO2genes{ $goserial } ){
            $genes = join(",", sort keys %{  $hGO2genes{ $goserial }{ $opts{t} } });
        }
        $goserial = "GO:" . $goserial;
        print OUT join("\t", $goserial, @arr, $genes), "\n";
    }
}
close IN;
close OUT;