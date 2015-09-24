#author         :Weihua Chen
#email          :chenwh550@hotmail.com
#last modified  :2006-1-9

#-------------------------------------------------------------------------
#please consult DOCUMENTATION for detailed information...
use strict;
use diagnostics;
use warnings;
use Getopt::Long;

my %opts=();
GetOptions(\%opts,"i:s","o:s","d:s");

if (!$opts{i} or !$opts{o} or !$opts{d}){
    print "----------------------------------------------------------------------
    USAGE: perl $0
        -i input eblast file
        -o out file
	-d TGI unigene2go database
----------------------------------------------------------------------\n";
    exit;
    
}
my $ifeblastn=$opts{i};
my $ofgo=$opts{o};
my $database=$opts{d};

open (EBLAST,$ifeblastn)||die "Cannot open file : $ifeblastn!!\n";
open (OUTGO,">$ofgo")||die "Cannot create output file : $ofgo!!\n" ;
open (TGIGO,$database)||die "Cannot open database : $database!!\n";
my %storehash;
my $id='';
while(<TGIGO>){
    if (/>(\S+)/){
	$id=$1;
    }elsif(/^(GO\:\d+)/){
	$storehash{$id}{$1}=1;
    }
}
close TGIGO;

my %hash;
while(<EBLAST>){
    chomp;
    next if (!/\d+/);
    next if (!$_);
    my ($q_name,$letter,$queryX,$queryY,$sbjctX,$sbjctY,$length,$score,$e_value,$overlap_total,$identity,$sbject,$anontation)=split(/\t/,$_);
    next if(exists $hash{$q_name});#this ensures every unigene counted only once
    if(exists $storehash{$sbject})
    {
	$hash{$q_name}=1;
	print OUTGO $q_name,"\t",join("\t",keys %{$storehash{$sbject}}),"\n";
    }
}
close EBLAST;
close OUTGO;

